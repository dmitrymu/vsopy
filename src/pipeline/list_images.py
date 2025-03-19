import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
import astropy.units as u
import ccdproc as ccdp
import numpy as np

from astropy.table import QTable
from astropy.time import Time
from vso import util

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from images for the session'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-i', '--image-dir', type=str, required=True, default=None, help='Image directory')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    img_layout = util.ImageLayout(args.image_dir)

    KEYS = set(['file', 'exptime', 'ccd-temp', 'filter', 'airmass', 'date-obs'])

    files = QTable([dict({key: row[key]
                          for key in row.colnames
                          if key in KEYS},
                         dir=d)
                    for d in img_layout.get_images(session).lights_dir.iterdir()
                    if d.is_dir()
                    for row in ccdp.ImageFileCollection(d).summary])
    files['image_id'] = [n+1 for n in range(len(files))]
    files['time'] = Time(files['date-obs'])
    files['exposure'] = files['exptime'] * u.second
    files['temperature'] = files['ccd-temp'] * u.deg_C
    files.add_column([str(row['dir'] / row['file']) for row in files], name='path')

    images = files['image_id', 'filter', 'time', 'exposure', 'airmass', 'temperature', 'path']
    images.meta = {'object': args.object, 'start': np.min(images['time']), 'finish': np.max(images['time'])}

    images.write(work_layout.get_session(session).images_file_path,
                format='ascii.ecsv', overwrite=args.overwrite)

    return 0

# Example: python3 list_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())