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
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--object', type=str, required=True, help='Object to be worked')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-c', '--common_dir', type=str, required=True, help='File tree root')

    return parser.parse_args()

def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object

    session_dir = Path(args.common_dir) /  Path('session') / object_dir
    if not session_dir.exists():
        session_dir.mkdir(parents=True)
    charts_dir = Path(args.common_dir) / 'charts_'

    data = vso.data.StarData(charts_dir)
    name =  args.object.replace('_', ' ')
    data.get_chart(name, 60*u.arcmin, maglimit=16*u.mag).write(session_dir / 'chart.ecsv', format='ascii.ecsv')

    return 0

# Example: python3 get_chart.py -O RR_Lyr -t 20230704 -c /srv/public

if __name__ == '__main__':
    sys.exit(main())