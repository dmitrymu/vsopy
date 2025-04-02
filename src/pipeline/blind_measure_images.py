import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
import astropy.units as u
import concurrent.futures as cf
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.nddata import CCDData
from astropy.table import QTable, Column, join, vstack
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from vso import phot
from vso import reduce
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from session images'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-p', '--parallel', type=int, default=1, help='Number of parallel workers')
    parser.add_argument('--snr', type=float, default=15, help='SNR threshold in dB')
    parser.add_argument('--max-separation', type=float, default=1.4, help='Star separation tolerance')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite output files')

    return parser.parse_args()

CONTEXT = None

def make_context(args):
    print('making context')
    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)
    solver = lambda path: reduce.astap_solver(path, session_layout.solved_dir)
    matcher = reduce.CalibrationMatcher(work_layout.calibr_dir)
    settings = util.Settings(session_layout.settings_file_path)
    global CONTEXT
    CONTEXT = (matcher, solver, settings.aperture)

def find_image_centroids(image, fwhm=10., threshold=5.):
    _, _, std = sigma_clipped_stats(image.data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
    sources = daofind(image.data)
    return QTable(dict(radec2000=image.wcs.pixel_to_world(sources['xcentroid'],
                                                          sources['ycentroid']),
                       auid=sources['id']))

def blind_measure_image(id, path, snr_th):
    print(f'measure {path}')
    global CONTEXT
    matcher, solver, aperture = CONTEXT
    result = phot.process_image(path, matcher, solver, find_image_centroids, aperture)
    result = result[result['snr'].value > snr_th]
    result['image_id'] = id
    return result['image_id', 'auid', 'radec2000', 'M', 'flux', 'snr', 'peak']

def build_star_table(tables, dist_th):
    stars = None

    id_map = QTable(dict(
        star_id = Column([], dtype='int32'),
        image_id = Column([], dtype='int32'),
        auid = Column([], dtype='int32')

    ))

    last_sid = 1

    for table in tables:
        c = stars['radec2000'] if stars is not None else SkyCoord(ra=[], dec=[], unit=(u.deg, u.deg))
        # look for matching stars (fast with K-D tree)
        x1, xs, sep, _ = table['radec2000'].search_around_sky(c, dist_th)

        # add new stars to the table
        mask = np.ones(len(table), dtype=bool)
        mask[xs]=False
        new_rows = table[mask]
        new_size = len(new_rows)
        star_ids = range(last_sid, last_sid + new_size)
        id_map = vstack([id_map, QTable(dict(
            star_id = star_ids,
            image_id = new_rows['image_id'],
            auid = new_rows['auid'],
        ))])
        new = QTable(dict(
            star_id = star_ids,
            radec2000=new_rows['radec2000'],
            err=np.zeros(new_size, dtype='float32')*u.arcsec,
            count= np.ones(new_size, dtype='int32')
        ))
        last_sid += new_size
        if (len(new) > 0):
            stars = vstack([stars, new]) if stars is not None else new

        # modify existing stars in the table
        if (len(xs) > 0):
            mask = np.ones(len(stars), dtype=bool)
            mask[x1]=False
            update = table[xs]
            old = stars[x1]
            id_map = vstack([id_map, QTable(dict(
                star_id = old['star_id'],
                image_id = update['image_id'],
                auid = update['auid'],
            ))])
            dra, ddec = old['radec2000'].spherical_offsets_to(update['radec2000'])
            weight = 1.0/(old['count']+1)
            old['radec2000'] = old['radec2000'].spherical_offsets_by(weight*dra, weight*ddec)
            old['err'] = np.sqrt(old['err']**2/old['count'] + sep**2)
            old['count'] += 1
            stars = vstack([stars[mask], old])
    mask = stars['err'] == 0.
    update = stars[mask]
    update['err'] = np.sqrt(dist_th**2/2)
    stars = vstack([stars[~mask], update])
    return stars, id_map

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)

    images = QTable.read(session_layout.images_file_path)

    blacklist = util.Blacklist(session_layout.blacklist_file_path)

    with cf.ProcessPoolExecutor(initializer=make_context,
                                initargs=(args,),
                                max_workers=args.parallel) as executor:
        futures = [(image['path'], executor.submit(blind_measure_image, image['image_id'], image['path'], args.snr))
                   for image in images if not blacklist.contains(image['path'])]

    def get_result(image_result):
        path = None
        try:
            path, future = image_result
            return future.result()
        except Exception as e:
            print(e)
            blacklist.add(path, e)
            return None

    tables = list([t for t in [get_result(f) for f in futures] if t is not None])
    stars, id_map = build_star_table(tables, args.max_separation*u.arcsec)
    result = join(vstack(tables), id_map, ['image_id', 'auid'])
    result.remove_column('auid')

    result.write(session_layout.root_dir / 'blind_measured.ecsv', format='ascii.ecsv', overwrite=args.overwrite)
    stars.write(session_layout.root_dir / 'detected_stars.ecsv', format='ascii.ecsv', overwrite=args.overwrite)
    blacklist.save(session_layout.blacklist_file_path)

    return 0

# Example: python3 measure_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())