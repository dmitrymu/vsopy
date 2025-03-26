import astropy.units as u
import numpy as np

from .. import phot
from ..util import MagErr, ValErr, MagErrDtype, ValErrDtype
from astropy.table import QTable, Column, join
from collections import namedtuple
from itertools import combinations
from scipy import linalg as sla

ClassicDiffTransform = namedtuple('ClassicDiffTransform', ['band', 'Ta', 'ka', 'Tb', 'kb'])

def build_classic_regression_input(session_layout, bands):
    """_summary_

    Args:
        name (_type_): _description_
        tag (_type_): _description_
        bands (_type_): _description_

    Returns:
        _type_: _description_
    """
    # A = [C1-C2, -X*(C1-C2)] = [dC, -X*dC];
    #  x = [T, k2];
    #  B = [(M1-M2) - (m1-m2)] = [dM - dm]

    def build_diff_column(data, col1, col2):
        coldata = [((row[col1]['mag'] - row[col2]['mag']).value,
                    np.sqrt(row[col1]['err']**2 + row[col2]['err']**2).value)
                    for row in data]
        return Column(coldata, unit=u.mag, dtype=MagErrDtype)

    def build_diff_sequence(sequence, selector, bands):
        b01, b02 = f'{bands[0]}_1', f'{bands[0]}_2'
        b11, b12 = f'{bands[1]}_1', f'{bands[1]}_2'
        A11 = QTable({
            'auid_1': sequence['auid'],
            b01: sequence[bands[0]],
            'C_1': sequence['C']
        })
        A12 = QTable({
            'auid_2': sequence['auid'],
            b02: sequence[bands[0]],
            'C_2': sequence['C']
        })
        A21 = QTable({
            'auid_1': sequence['auid'],
            b11: sequence[bands[1]]
        })
        A22 = QTable({
            'auid_2': sequence['auid'],
            b12: sequence[bands[1]]
        })
        A0 = join(join(join(join(selector, A11, 'auid_1'), A12, 'auid_2'), A21, 'auid_1'), A22, 'auid_2')
        A0[f'd{bands[0]}'] = build_diff_column(A0, b01, b02)
        A0[f'd{bands[1]}'] = build_diff_column(A0, b11, b12)
        A0['dC'] = build_diff_column(A0, 'C_1', 'C_2')
        return A0['auid_1', 'auid_2', f'd{bands[0]}', f'd{bands[1]}', 'dC']

    def build_diff_measurement(measured, selector, bands):
        b01, b02 = f'{bands[0]}_1', f'{bands[0]}_2'
        b11, b12 = f'{bands[1]}_1', f'{bands[1]}_2'
        m11 = QTable({
            'batch_id': measured['batch_id'],
            'auid_1': measured['auid'],
            b01: measured[f'instr {bands[0]}']
        })
        m12 = QTable({
            'batch_id': measured['batch_id'],
            'auid_2': measured['auid'],
            b02: measured[f'instr {bands[0]}']
        })
        m21 = QTable({
            'batch_id': measured['batch_id'],
            'auid_1': measured['auid'],
            b11: measured[f'instr {bands[1]}']
        })
        m22 = QTable({
            'batch_id': measured['batch_id'],
            'auid_2': measured['auid'],
            b12: measured[f'instr {bands[1]}']
        })

        m = join(join(join(join(selector, m11, ['batch_id', 'auid_1']), m12, [
                 'batch_id', 'auid_2']), m21, ['batch_id', 'auid_1']), m22, ['batch_id', 'auid_2'])
        m[f'instr d{bands[0]}'] = build_diff_column(m, b01, b02)
        m[f'instr d{bands[1]}'] = build_diff_column(m, b11, b12)
        return m['batch_id', 'auid_1', 'auid_2', f'instr d{bands[0]}', f'instr d{bands[1]}']

    provider = phot.BatchDataProvider(session_layout)

    sequence = provider.sequence_band_pair(bands)
    sequence['C'] = build_diff_column(sequence, bands[0], bands[1]) #Column(sequence[bands[0]]['mag'] - sequence[bands[1]]['mag'])

    id1, id2 = zip(*(combinations(sequence['auid'], 2)))
    # stars = QTable([sequence['auid']])
    # cartesian = join(stars, stars, join_type='cartesian')
    cartesian = QTable(dict(auid_1 = id1, auid_2 = id2))
    sequence_selector = cartesian[cartesian['auid_1'] != cartesian['auid_2']]
    A0 = build_diff_sequence(sequence, sequence_selector, bands)

    batches = QTable(dict(batch_id = provider.batches_['batch_id']))
    batch_sequence_selector = join(batches, sequence_selector, join_type='cartesian')
    measured = provider.batch_band_pair(bands)
    diff_measured = build_diff_measurement(measured, batch_sequence_selector, bands)

    airmass = QTable(dict(
        batch_id = provider.batches_['batch_id'],
        airmass = Column([(row['airmass'], row['airmass_range']) for row in provider.batches_], dtype=ValErrDtype)
        ))

    result = join(join(diff_measured, A0, ['auid_1', 'auid_2']), airmass, 'batch_id')
    return result


def calc_classic_diff_transform(data, bands):
    def calc_band_transform(data, band, A, Qjj):
        B = data[f'd{band}']['mag'] - data[f'instr d{band}']['mag']
        X, norm2, _, _ = sla.lstsq(A, B) # [Ta, ka], ...
        # see https://en.wikipedia.org/wiki/Ordinary_least_squares#Intervals
        err = np.sqrt( Qjj * norm2 / (len(data)-2)) # [Ta_err, ka_err]
        return X, err

    A = np.array([data['dC']['mag'], -data['airmass']['val'] * data['dC']['mag']]).T
    Qjj = np.diag(np.linalg.inv( np.matmul(A.T, A)))
    X1, err1 = calc_band_transform(data, bands[0], A, Qjj)
    X2, err2 = calc_band_transform(data, bands[1], A, Qjj)
    return ClassicDiffTransform(bands,
                                ValErr(X1[0], err1[0]), ValErr(X1[1], err1[1]),
                                ValErr(X2[0], err2[0]), ValErr(X2[1], err2[1]))

def calc_classic_diff_transform_weighted(data, bands):
   def calc_band_transform_weighted(data, band, A, errA):
      errB = np.sqrt(data[f'd{band}']['err']**2 + data[f'instr d{band}']['err']**2 + errA)
      # see https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix
      W = np.diag(np.reciprocal(errB))
      Aw = np.dot(W, A)
      Qjj = np.diag(np.linalg.inv(np.matmul(Aw.T, Aw)))
      B = data[f'd{band}']['mag'] - data[f'instr d{band}']['mag']
      Bw = np.dot(B, W)
      X, res, _, _ = sla.lstsq(Aw, Bw) # [Ta, ka], ...
      err = np.sqrt( Qjj * res / (len(data)-2)) # [Ta_err, ka_err]
      return X, err

   A = np.array([data['dC']['mag'], -data['airmass']['val'] * data['dC']['mag']]).T
   errA = data[f'dC']['err']**2 + (data[f'dC']['mag'] * data['airmass']['val'])**2 \
                                 * ((data[f'dC']['err']/data[f'dC']['mag'])**2 \
                                    + (data['airmass']['err']/data['airmass']['val'])**2)
   X1, err1 = calc_band_transform_weighted(data, bands[0], A, errA)
   X2, err2 = calc_band_transform_weighted(data, bands[1], A, errA)
   return ClassicDiffTransform(bands,
                                ValErr(X1[0], err1[0]), ValErr(X1[1], err1[1]),
                                ValErr(X2[0], err2[0]), ValErr(X2[1], err2[1]))
