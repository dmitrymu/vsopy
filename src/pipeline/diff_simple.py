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
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--object', type=str,
                        required=True, help='Object to be worked')
    parser.add_argument('-t', '--tag', type=str,
                        required=True, help='Tag (date)')
    parser.add_argument('-c', '--common_dir', type=str,
                        required=True, help='File tree root')
    parser.add_argument('--observer', type=str,
                        required=True, help='AAVSO observer code')

    return parser.parse_args()


def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object  # '20240725/SA38'

    session_dir = args.common_dir / Path('session') / object_dir

    with open(session_dir / 'settings.json') as file:
       settings = json.load(file)

    provider = phot.DataProvider(
        QTable.read(session_dir / 'photometry.ecsv'),
        QTable.read(session_dir / 'chart.ecsv')
    )

    def total_err(table, band):
        return np.sqrt(np.sum(table[band]['err']**2)/len(table))

    xfm = phot.BatchTransformer(('B', 'V'), settings["diff_photometry"])
    bv = xfm.calculate(provider)
    V_BV_err = total_err(bv, 'V')
#    bv.write(session_dir / 'diff_simple_bv.ecsv', overwrite=True)

    xfm = phot.BatchTransformer(('V', 'Rc'), settings["diff_photometry"])
    vr = xfm.calculate(provider)
    V_VR_err = total_err(vr, 'V')
    R_VR_err = total_err(vr, 'Rc')
#     result.write(session_dir / 'verify_simple_vr.ecsv', overwrite=True)

    xfm = phot.BatchTransformer(('Rc', 'Ic'), settings["diff_photometry"])
    ri = xfm.calculate(provider)
    R_RI_err = total_err(ri, 'Rc')
#     result.write(session_dir / 'verify_simple_ri.ecsv', overwrite=True)

    with open(session_dir/'report-simple.txt', mode='w') as f:
       report = data.AavsoReport(f,
                                 provider.target_name,
                                 provider.chart_id,
                                 args.observer)
       report.header()
       report.body(bv, 'B')
       report.body(bv if V_BV_err < V_VR_err else vr, 'V')
       report.body(vr if R_VR_err < R_RI_err else vr, 'Rc')
       report.body(ri, 'Ic')

    return 0

# Example: python3 diff_simple.py -O RR_Lyr -t 20230704 -c /srv/public --observer XYZ


if __name__ == '__main__':
    sys.exit(main())
