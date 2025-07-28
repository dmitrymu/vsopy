import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import numpy as np
import argparse
import json
from vsopy import phot, data, util
from pathlib import Path
from astropy.table import QTable
from collections import namedtuple
from operator import attrgetter

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

    session = util.Session(tag=args.tag, name=args.object)
    layout = util.WorkLayout(args.work_dir)
    session_layout = layout.get_session(session)
    settings = util.Settings(session_layout.settings_file_path)

    provider = phot.DataProvider(
        QTable.read(session_layout.photometry_file_path),
        QTable.read(session_layout.chart_file_path)
    )

    def reject_outliers(data, band, sigma=2.0):
        dc = data[f'check {band}']['mag'] - data[f'check std {band}']['mag']
        dc_mean = np.mean(dc)
        dc_std = np.std(dc)
        flt = np.abs(dc - dc_mean) < sigma*dc_std
        return data[flt]

    def total_err(table, band):
        return np.sqrt(np.sum(table[band]['err']**2)/len(table))

    BandResult = namedtuple('BandResult', ['err', 'data'])

    bands = util.band_pairs(settings.bands)
    result = {}
    for band in bands:
        xfm = phot.BatchTransformer(
            band,
            settings.get_comp(band),
            settings.get_check(band))

        ab = xfm.calculate(provider)
        a_err = total_err(ab, band[0])
        b_err = total_err(ab, band[1])
        result.setdefault(band[0], []).append(BandResult(err=a_err, data=ab))
        result.setdefault(band[1], []).append(BandResult(err=b_err, data=ab))

    with open(session_layout.root_dir / 'report-simple.txt', mode='w') as f:
      report = data.AavsoReport(f,
                                provider.target_name,
                                provider.chart_id,
                                args.observer)
      report.header()
      for band in settings.bands:
        dt = min(result[band], key=attrgetter('err')).data
        report.body(reject_outliers(dt, band), band)

    return 0

# Example: python3 diff_simple.py -O RR_Lyr -t 20230704 -w /home/user/work --observer XYZ


if __name__ == '__main__':
    sys.exit(main())
