import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import argparse
from vso import phot
from pathlib import Path
from astropy.table import QTable, Column, vstack

import astropy.units as u

from itertools import chain

class Transformer:
    def __init__(self, band) -> None:
        self.band_ = band

    def calculate(self, provider):
        batches = provider.get_batches(self.band_)

        result = [(batch,
                  phot.SimpleTransform.create(batch, self.band_[0], self.band_[1]))
                  for batch in batches]
        return result

    def validate(self, provider):

        def validate_batch(batch, xfm):

            N = len(batch)

            def split(batch, n):
                comp = batch[n]
                target_filter = [True]*N
                target_filter[n] = False
                targets = batch[target_filter]
                return comp, targets


            batch_trials = [split(batch, n) for n in range(N)]
            trials = chain.from_iterable([[(comp, target, xfm(target, comp))
                     for target in targets]
                    for comp, targets in batch_trials])

            columns = ['auid', self.band_[0], self.band_[1],
                           f"instr {self.band_[0]}", f"instr {self.band_[1]}"]
            def mags_to_row(x):
                t = QTable([
                    Column([(x[0][0].value, x[0][1].value)],
                           dtype=[('mag', 'f4'), ('err', 'f4')],
                           unit=u.mag,
                           name='A'),
                    Column([(x[1][0].value, x[1][1].value)],
                           dtype=[('mag', 'f4'), ('err', 'f4')],
                           unit=u.mag,
                           name='B')
                ])

                return list(t[0]['A', 'B'])
            cA = f"pred  {self.band_[0]}"
            cB = f"pred  {self.band_[1]}"
            tt = QTable(rows=[list(t[0][columns])
                              + list(t[1][columns])
                              + mags_to_row(t[2])
                              for t in trials],
                        names=[f"comp {f}" for f in columns]
                        + [f"targ {f}" for f in columns]
                        + [cA, cB])
            return tt
        return vstack([validate_batch(batch, xfm) for batch, xfm in self.calculate(provider)])



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--object', type=str, required=True, help='Object to be worked')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-c', '--common_dir', type=str, required=True, help='File tree root')

    return parser.parse_args()

def main():
    args = parse_args()

    object_dir = Path(args.tag) / args.object #'20240725/SA38'

    session_dir = args.common_dir /  Path('session') / object_dir
    provider = phot.DataProvider(
        QTable.read(session_dir / 'photometry.ecsv'),
        QTable.read(session_dir / 'chart.ecsv')
    )

    xfm = Transformer(('B', 'V'))
    result = xfm.validate(provider)
    result.write(session_dir / 'verify_simple_bv.ecsv', overwrite=True)

    xfm = Transformer(('V', 'Rc'))
    result = xfm.validate(provider)
    result.write(session_dir / 'verify_simple_vr.ecsv', overwrite=True)

    xfm = Transformer(('Rc', 'Ic'))
    result = xfm.validate(provider)
    result.write(session_dir / 'verify_simple_ri.ecsv', overwrite=True)

    return 0

# Example: python3 verify_simple_xfm.py -O RR_Lyr -t 20230704 -c /srv/public

if __name__ == '__main__':
    sys.exit(main())