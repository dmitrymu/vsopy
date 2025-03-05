import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
import vso.data
import vso.util
from astropy.table import unique, vstack

import astropy.units as u
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description='download photometric sequence to the session directory')
    parser.add_argument('-O', '--object', type=str,
                        required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str,
                        required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str,
                        required=True, help='Work directory')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = vso.util.Session(tag=args.tag, name=args.object)
    layout = vso.util.WorkLayout(args.work_dir)
    session_layout = layout.get_session(session)

    data = vso.data.StarData(layout.charts_dir)
    name =  args.object.replace('_', ' ')
    if data.is_std_field(name):
        objects = [n for n in data.std_fields['name']
                   if n == name or n.startswith(f"{name} ")]

        stars = unique(vstack([data.get_chart(name)
                               for name in objects]))
        stars.write(session_layout.chart_file_path,
                    format='ascii.ecsv',
                    overwrite=args.overwrite)
    else:
        data.get_chart(name,
                       60*u.arcmin,
                       maglimit=16 * u.mag).write(session_layout.chart_file_path,
                                                  format='ascii.ecsv',
                                                  overwrite=args.overwrite)

    return 0

# Example: python3 get_chart.py -O RR_Lyr -t 20230704 -w /home/user/work

if __name__ == '__main__':
    sys.exit(main())