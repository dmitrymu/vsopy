import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
import astropy.units as u
import ccdproc as ccdp

from astropy.table import QTable
import numpy as np
from vso import util
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract instrumental magnitudes from images for the session'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    # parser.add_argument('-i', '--image-dir', type=str, required=True, default=None, help='Image directory')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    session_layout = util.WorkLayout(args.work_dir).get_session(session)

    images = QTable.read(session_layout.images_file_path)
    images.sort('time')
    bands = set(images['filter'])
    batch_size = len(bands)
    order = list(images[0:4]['filter'])

    images.add_column(0, name='batch_id')
    images.add_column(0*u.second, name='time_range')
    images.add_column(0*images['temperature'].unit, name='temperature_range')
    images.add_column(0., name='airmass_range')
    batches = images['batch_id', 'time', 'temperature', 'airmass',
                     'time_range', 'temperature_range', 'airmass_range'][:0].copy()
    batch_images = images['batch_id', 'image_id', 'filter'][:0].copy()

    next = 0
    id = 0
    while next < len(images) - batch_size + 1:
        batch = images[next:next+batch_size]
        if np.all(batch['filter'] == order):
            id += 1
            d = np.max(batch['airmass']) - np.min(batch['airmass'])
            batches.add_row({
                'batch_id': id,
                'time': np.mean(batch['time']),
                'time_range': np.max(batch['time']) - np.min(batch['time']),
                'temperature': np.mean(batch['temperature']),
                'temperature_range': np.max(batch['temperature']) - np.min(batch['temperature']),
                'airmass': np.mean(batch['airmass']),
                'airmass_range': d #np.max(batch['airmass']) - np.min(batch['airmass'])
            })
            for image in batch:
               batch_images.add_row(dict(
                   batch_id=id,
                   image_id=image['image_id'],
                   filter=image['filter']
               ))
            next += batch_size
        else:
            print(f"Skipped file #{next+1} '{images['file']}'")
            next += 1

    batches.meta = {'object': args.object, 'start': np.min(images['time']), 'finish': np.max(images['time'])}

    batches.write(session_layout.root_dir / 'batches.ecsv',
                  format='ascii.ecsv', overwrite=args.overwrite)

    batch_images.write(session_layout.root_dir / 'batch_images.ecsv',
                       format='ascii.ecsv', overwrite=args.overwrite)

    return 0

# Example: python3 batch_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())