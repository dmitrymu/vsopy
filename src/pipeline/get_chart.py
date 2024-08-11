import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import sys
import argparse
import vso.data
from pathlib import Path
from astropy.table import QTable

import astropy.units as u
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description='download photometric sequence to the session directory')
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')

    return parser.parse_args()

def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object

    session_dir = Path(args.work_dir) /  Path('session') / object_dir
    if not session_dir.exists():
        session_dir.mkdir(parents=True)
    charts_dir = Path(args.work_dir) / 'charts_'

    data = vso.data.StarData(charts_dir)
    name =  args.object.replace('_', ' ')
    data.get_chart(name, 60*u.arcmin, maglimit=16*u.mag).write(session_dir / 'chart.ecsv', format='ascii.ecsv')

    return 0

# Example: python3 get_chart.py -O RR_Lyr -t 20230704 -w /home/user/work

if __name__ == '__main__':
    sys.exit(main())