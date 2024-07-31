import astropy.units as u
import ccdproc as ccdp
import numpy as np
from .. import util
from ..data import CameraRegistry
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable, Column
from astropy.wcs import WCS
from photutils.aperture import (ApertureStats,
                                CircularAperture,
                                SkyCircularAperture,
                                SkyCircularAnnulus)


def measure_photometry(image, stars, r, r_in=None, r_out=None):
    result = QTable(stars['auid', 'sky_centroid'])
    if 'name' in stars.colnames:
        result['name'] = stars['name']

    _, median, _ = sigma_clipped_stats(image.data, sigma=3.0)
    corrected = CCDData(image)
    corrected.data = corrected.data - median

    apr = SkyCircularAperture(stars['sky_centroid'], r=r)
    ap_stats = ApertureStats(corrected, apr)
    ann = SkyCircularAnnulus(stars['sky_centroid'],
                             r_in=r_in if r_in is not None else r * 2,
                             r_out=r_out if r_out is not None else r*3)
    ann_stats = ApertureStats(corrected, ann)

    zero_level = 1_000_000 * u.electron / u.second
    exp = image.header['EXPTIME'] * u.second
    Ft = ap_stats.sum / exp
    Fb = ap_stats.sum_aper_area.value * ann_stats.mean / exp
    FtFb = Ft - Fb
    Ft2Fb = Ft + 2 * Fb
    result['flux'] = FtFb
    result['snr'] = 10 * np.log10(FtFb.value/np.sqrt(Ft2Fb.value)) * u.db
    mag = -2.5 * np.log10(result['flux'] / zero_level).value

    Ft_err = ap_stats.sum_err
    Fb_err = ap_stats.sum_aper_area.value * (ann_stats.sum_err / ann_stats.sum_aper_area.value)
    flux_err = np.sqrt(Ft_err*Ft_err + Fb_err*Fb_err) / exp
    mag_err = (2.5 * flux_err / result['flux'] / np.log(10)).value
    result['M'] = Column(list(zip(mag, mag_err)),
                         unit=u.mag,
                         dtype=[('mag','f4'), ('err', 'f4')])

    camera = CameraRegistry.get(image.header['instrume'])
    result['peak'] = np.nan if not camera else ap_stats.max.value / camera.max_adu.value

    return result
