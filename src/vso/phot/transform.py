import astropy.units as u
import numpy as np
from scipy import stats as sst
from collections import namedtuple
from astropy.table import QTable, Column, join

ValErr = namedtuple('ValErr', ['val', 'err'])
ValErrDtype = [('val', 'f4'), ('err', 'f4')]
MagErr = namedtuple('ValErr', ['mag', 'err'])
MagErrDtype = [('mag', 'f4'), ('err', 'f4')]
SimpleTransform = namedtuple('SimpleTransform', ['Ta', 'Tb', 'Tab'])

def create_simple_transform(A, B, a, b):
    AB = A - B
    ab = a - b
    Aa =A - a
    Bb = B - b
    reg_ab = sst.linregress(AB, ab)
    reg_Aa = sst.linregress(AB, Aa)
    reg_Bb = sst.linregress(AB, Bb)

    if reg_Aa.stderr == 0 or reg_Bb.stderr == 0 or reg_ab.stderr == 0:
        raise Exception("Zero stderr")

    result = SimpleTransform(
        ValErr(reg_Aa.slope, reg_Aa.stderr),
        ValErr(reg_Bb.slope, reg_Bb.stderr),
        ValErr(1/reg_ab.slope, reg_ab.stderr))
    return result

def apply_simple_transform(xfm, A_c, B_c, a_c, b_c, a_t, b_t):
    def transform(T_a, T_ab, A_c, B_c, a_c, b_c, a_t, b_t):
        Ta, Ta_err = T_a
        Tab, Tab_err = T_ab
        Ac, Ac_err = A_c
        Bc, Bc_err = B_c
        ac, ac_err = a_c
        bc, bc_err = b_c
        at, at_err = a_t
        bt, bt_err = b_t

        # c_t = a_t - b_t
        atbt = at-bt
        # c_C = a_c - b_c
        acbc = ac-bc

        atbtacbc = atbt - acbc
        atbtacbc_err = np.sqrt(at_err**2 + bt_err**2 + ac_err**2 + bc_err**2)
        Tab_abab = Tab * (atbt - acbc)
        Tab_abab_err = Tab * atbtacbc * \
            np.sqrt((Tab_err/Tab)**2 + (atbtacbc_err/atbtacbc)**2)
        # (1)
        AtBt = (Ac-Bc) + Tab * (atbt - acbc)
        AtBt_err = np.sqrt(Ac_err**2 + Bc_err**2 + Tab_abab_err**2)
        # (2)
        Acac = Ac - ac
        Ta_Tab_err = Ta * Tab_abab * \
            np.sqrt((Ta_err/Ta)**2 + (Tab_abab_err/Tab_abab)**2)
        At = at + Acac + Ta * Tab * (atbt - acbc)
        At_err = np.sqrt(at_err**2 + Ac_err**2 + ac_err**2 + Ta_Tab_err**2)

        # (3)
        Bt = At - AtBt
        Bt_err = np.sqrt(At_err**2 + AtBt_err**2)

        return MagErr(At, At_err), MagErr(Bt, Bt_err)

    Ta, Tb, Tab = xfm
    r1 =  transform(Ta, Tab, A_c, B_c, a_c, b_c, a_t, b_t)
    r2 =  transform(Tb, Tab, B_c, A_c, b_c, a_c, b_t, a_t)
    return r1[0], r2[0], r2[1], r1[1]

def batch_create_simple_transform(provider, bands):

    data = provider.batch_and_sequence_band_pair(bands)
    grouped = data.group_by(['batch_id'])

    xfm = [create_simple_transform(batch[bands[0]]['mag'], batch[bands[1]]['mag'],
                                   batch[provider.instr(bands[0])]['mag'], batch[provider.instr(bands[1])]['mag'])
            for batch in grouped.groups]
    xfm_dtype = [('Ta', ValErrDtype), ('Tb', ValErrDtype), ('Tab', ValErrDtype)]
    return QTable({
        'batch_id': grouped.groups.keys['batch_id'],
        'xfm': Column(xfm, dtype = xfm_dtype)
    })

def batch_diff_photometry(provider, bands, comparison_auid, target_auid=None):


    xfm_table = batch_create_simple_transform(provider, bands)

    comp = provider.batch_comp_star(bands, comparison_auid)
    targ = provider.batch_target_star(bands, target_auid if target_auid is not None else provider.target_auid)
    xfm_input = join(join(comp, targ, 'batch_id'), xfm_table, 'batch_id')
    transformed = [apply_simple_transform(row['xfm'], row[bands[0]].value, row[bands[1]].value,
                            row[provider.instr(bands[0])].value, row[provider.instr(bands[1])].value,
                            row[provider.targ(bands[0])].value, row[provider.targ(bands[1])].value
                            ) for row in xfm_input]
    a, b, _, _ =zip(*transformed)
    return QTable({
        'batch_id': xfm_input['batch_id'],
        bands[0]: Column(list(a), dtype = MagErrDtype, unit=u.mag),
        bands[1]: Column(list(b), dtype = MagErrDtype, unit=u.mag)
        })
