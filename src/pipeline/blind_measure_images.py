import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
import astropy.units as u
import concurrent.futures as cf
import quads

from astropy.coordinates import SkyCoord
from astropy.nddata import CCDData
from astropy.table import QTable, Column, vstack
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

def blind_measure_image(id, path):
    print(f'measure {path}')
    global CONTEXT
    matcher, solver, aperture = CONTEXT

    image = reduce.update_wcs(CCDData.read(path, unit='adu'), solver(path))
    calibration = matcher.match(image.header)
    reduced = reduce.calibrate_image(image,
                                    dark=calibration.dark,
                                    flat=calibration.flat)
    print('Searching stars')
    _, _, std = sigma_clipped_stats(reduced.data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=10.0, threshold=5.*std)
    sources = daofind(reduced.data)
    centroids = QTable(dict(radec2000=reduced.wcs.pixel_to_world(sources['xcentroid'], sources['ycentroid']),
                        auid=sources['id']))
    print(f'{len(sources)} found')

    result = phot.measure_photometry(reduced, centroids, aperture)
    result.remove_column('radec2000')
    result.rename_column('sky_centroid', 'radec2000')
    result['image_id'] = id
    return result['image_id', 'auid', 'radec2000', 'M', 'flux', 'snr', 'peak']

def aggregate_tables(tables, snr_th, dist_th):
    
    result = QTable(tables[0])[[]] # create from template
    result.add_column(Column([], dtype='int32'), name='star_id')
    result.remove_column('auid')
    #result.remove_column('radec2000')
    star_id = 0
    tree = quads.QuadTree((180.0, 0.0), 360.0, 180.0)    
    for t in tables:
        if t is not None:
            for star in t:
                #print(star)
                c = star['radec2000']
                if star['snr'].value > snr_th:
                    nn = tree.nearest_neighbors((c.ra.value, c.dec.value), count=1)
                    id = 0
                    if len(nn) < 1 or SkyCoord(ra=nn[0].x, dec=nn[0].y, unit=(u.deg, u.deg)).separation(c) > dist_th:
                        # print(f'creating star {star_id}, nearest {None if len(nn) < 1 else nn[0]}')
                        star_id += 1
                        tree.insert((c.ra.value, c.dec.value), data = star_id)
                        id = star_id
                    else:
                        id = nn[0].data
                    # print(f'adding star {star_id}')
                    result.add_row({
                        'image_id': star['image_id'],
                        'star_id': id,
                        'radec2000': star['radec2000'],
                        'M': star['M'],
                        'flux': star['flux'],
                        'snr': star['snr'],
                        'peak': star['peak']
                    })

    stars=QTable(dict(
                star_id =[0],
                radec2000=SkyCoord([0], [0], unit=(u.deg, u.deg))
                ))[[]]
    for node in tree:
        stars.add_row([
            node.data,
            SkyCoord(ra=node.x, dec=node.y, unit=(u.deg, u.deg))
        ])
    return result, stars



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
        futures = [(image['path'], executor.submit(blind_measure_image, image['image_id'], image['path']))
                  for image in images if not blacklist.contains(image['path'])]

    def get_result(image_result):
        path = None
        try:
            path, future = image_result
            return future.result()
        except Exception as e:
            blacklist.add(path, e)
            return None

    tables = list([get_result(f) for f in futures])
    tables = [t for t in tables if t is not None]
    if len(tables) < 1 or tables[0] is None:
        return 1

    #result = vstack([t for t in tables if t is not None])
    result, stars = aggregate_tables(tables, 15, .5*u.arcsec)

    result.write(session_layout.root_dir / 'blind_measured.ecsv', format='ascii.ecsv', overwrite=args.overwrite)
    stars.write(session_layout.root_dir / 'detected_stars.ecsv', format='ascii.ecsv', overwrite=args.overwrite)
    blacklist.save(session_layout.blacklist_file_path)

    return 0

# Example: python3 measure_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())