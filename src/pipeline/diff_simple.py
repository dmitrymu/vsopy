import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import numpy as np
import argparse
import json
from vso import phot, data
from pathlib import Path
from astropy.table import QTable


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate differential photometry report for the session')
    parser.add_argument('-O', '--object', type=str,
                        required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str,
                        required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str,
                        required=True, help='Work directory')
    parser.add_argument('--observer', type=str,
                        required=True, help='AAVSO observer code')

    return parser.parse_args()


def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object  # '20240725/SA38'

    session_dir = args.work_dir / Path('session') / object_dir

    with open(session_dir / 'settings.json') as file:
       settings = json.load(file)

    provider = phot.DataProvider(
        QTable.read(session_dir / 'photometry.ecsv'),
        QTable.read(session_dir / 'chart.ecsv')
    )


    def get_comp(band, settings):
        band_settings = settings["diff_photometry"][f"{band[0]}{band[1]}"]
        return band_settings['comp']

    def get_check(band, settings):
        band_settings = settings["diff_photometry"][f"{band[0]}{band[1]}"]
        return  band_settings['check'] if 'check' in band_settings else None


    def total_err(table, band):
        return np.sqrt(np.sum(table[band]['err']**2)/len(table))

    band = ('B', 'V')
    xfm = phot.BatchTransformer(band,
                                get_comp(band, settings),
                                get_check(band, settings))
    bv = xfm.calculate(provider)
    V_BV_err = total_err(bv, 'V')

    band = ('V', 'Rc')
    xfm = phot.BatchTransformer(band,
                                get_comp(band, settings),
                                get_check(band, settings))
    vr = xfm.calculate(provider)
    V_VR_err = total_err(vr, 'V')
    R_VR_err = total_err(vr, 'Rc')

    band = ('Rc', 'Ic')
    xfm = phot.BatchTransformer(band,
                                get_comp(band, settings),
                                get_check(band, settings))
    ri = xfm.calculate(provider)
    R_RI_err = total_err(ri, 'Rc')

    with open(session_dir/'report-simple.txt', mode='w') as f:
       report = data.AavsoReport(f,
                                 provider.target_name,
                                 provider.chart_id,
                                 args.observer)
       report.header()
       report.body(bv, 'B')
       report.body(bv if V_BV_err < V_VR_err else vr, 'V')
       report.body(vr if R_VR_err < R_RI_err else ri, 'Rc')
       report.body(ri, 'Ic')

    return 0

# Example: python3 diff_simple.py -O RR_Lyr -t 20230704 -w /home/user/work --observer XYZ


if __name__ == '__main__':
    sys.exit(main())
