import numpy as np
from scipy import stats as sst
from typing import Tuple

ValErr = Tuple[float, float]

def transform(T_a: ValErr, T_ab: ValErr, A_c: ValErr, B_c: ValErr, a_c: ValErr, b_c: ValErr, a_t: ValErr, b_t: ValErr):
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

        return (At, At_err), (Bt, Bt_err)


class SimpleTransform:
    """Magnitude transformation for differential photometry.

    It uses simplified approach which ignores second order extinction.
    See Gary B.L., CCD Transformation Equations for Use with Single-Image Photometry.
    This transform must be calculated for each pair of images separately and applied
    to the instrumental magnitudes extracted from those images.
    """

    def __init__(self, band_a: str, band_b: str, Ta: ValErr, Tb: ValErr, Tab: ValErr):
        self.Ta, self.Ta_err = Ta
        self.Tb, self.Tb_err = Tb
        self.Tab, self.Tab_err = Tab
        self.band_a = band_a
        self.band_b = band_b

    def __repr__(self) -> str:
        return f"SimpleTransform(T_{self.band_a} = {self.Ta:.3g} +/- {self.Ta_err:.3g}; T_{self.band_a}{self.band_b} = {self.Tab:.3g} +/- {self.Tab_err:.3g})"

    def __call__(self, target, comp):
        """Calculate standard magnitude of the target using its instrumental magnitude, comparison star, and transform.

        For color bands A and B defined on transform creation, the following magnitudes are known:
        * a_t and b_t - instrumental for target star;
        * a_c and b_c - instrumental for comparison star;
        * A_c and B_c - standard for comparison star.

        Standard magnitudes of the target A_t and B_t are calculated as follows:

        (1) C_t = A_t - B_t = (A_c - B_c) + T_ab * ((a_t-b_t) - (a_c-b_c))
        (2) A_t = a_t + (A_c-a_c) + T_b * (C_t - (A_c - B_c))
        (3) B_t = A_t - C_t

        All magnitudes and transform coefficients have their uncertainties that should
        be propagated through the equations above.

        Args:
            target (dict-like): target star (instrumental)
            comp (dict_like): comparison star (instrumental and standard)
        """
        def split_mag(v):
            return (v['mag'], v['err'])

        Ac, Ac_err = split_mag(comp[self.band_a])
        Bc, Bc_err = split_mag(comp[self.band_b])
        ac, ac_err = split_mag(comp[f'instr {self.band_a}'])
        bc, bc_err = split_mag(comp[f'instr {self.band_b}'])
        at, at_err = split_mag(target[f'instr {self.band_a}'])
        bt, bt_err = split_mag(target[f'instr {self.band_b}'])

        r1 =  transform((self.Ta, self.Ta_err),
                         (self.Tab, self.Tab_err),
                         (Ac, Ac_err),
                         (Bc, Bc_err),
                         (ac, ac_err),
                         (bc, bc_err),
                         (at, at_err),
                         (bt, bt_err))

        r2 =  transform((self.Tb, self.Tb_err),
                         (self.Tab, self.Tab_err),
                         (Bc, Bc_err),
                         (Ac, Ac_err),
                         (bc, bc_err),
                         (ac, ac_err),
                         (bt, bt_err),
                         (at, at_err))
        return r1[0], r1[1], r2[1], r2[0]

    @staticmethod
    def create(batch, band_a: str, band_b: str):
        """Create transform from star photometry for two images

        Given an ensemble of stars with both standard (A, B) and instrumental
        (a,b) magnitudes known for two bands, we fit two linear regressions:

        (1)  a - b = T_ab * (A - B) + C_ab
        (2)  A-a = T_a * (A - B) + C_a
        (3)  B-b = T_b * (A - B) + C_b

        T_ab determines transformation of instrumental color index to standard
        color index.  T_a corrects transformation from instrumental to standard
        magnitude with respect to color index. See comment on
        SimpleTransform.__call__ for details.

        Args:
            batch (table-like): instrumental and standard magnitudes for image pair
            band_a (str): band (filter) of the first image
            band_b (str): band (filter) of the second image

        Returns:
            SimpleTransform: transform to calculate target standard magnitude
            given comparison star.
        """
        def no_nan_mag(batch, band_a, band_b):
            filter_nan = np.all([~np.isnan(batch[band_a]['mag']),
                                ~np.isnan(batch[band_b]['mag'])], axis=0)
            return batch[filter_nan]

        b = no_nan_mag(batch, band_a, band_b)
        AB = b[band_a]['mag'] - b[band_b]['mag']
        ab = b[f'instr {band_a}']['mag'] - b[f'instr {band_b}']['mag']
        Aa = b[band_a]['mag'] - b[f'instr {band_a}']['mag']
        Bb = b[band_b]['mag'] - b[f'instr {band_b}']['mag']
        reg_ab = sst.linregress(AB, ab)
        reg_Aa = sst.linregress(AB, Aa)
        reg_Bb = sst.linregress(AB, Bb)

        if reg_Aa.stderr == 0 or reg_ab.stderr == 0:
            raise Exception("Zero stderr%")

        result = SimpleTransform(band_a, band_b,
                                 (reg_Aa.slope, reg_Aa.stderr),
                                 (reg_Bb.slope, reg_Bb.stderr),
                                 (1/reg_ab.slope, reg_ab.stderr))
        return result
