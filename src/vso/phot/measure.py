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


def measure_photometry(image, stars, r_ap, r_ann):
    """Extract aperture photometry data from the image.

    Given the calibrated image and the list of star centroids in sky coordinates,
    measure star flux using circular aperture and sky annulus.  For detailed discussion
    see Warner B.D., A Practical Guide to Lightcurve Photometry and Analysis, 4.2.

    Flux counting is performed by photutils package
    (https://photutils.readthedocs.io/en/stable/aperture.html). Uncertainty from
    the image is propagated to the flux values. Flux is normalized by the image
    exposure to electrons per second.

    The flux for central aperture F_c contains both star and sky electrons, the flux
    for annulus F_a - sky electrons only. Number of sky electrons in the central
    aperture is.

    F_{sky} = F_{sky mean} * (pixel count for central aperture)
    F_{star} = F_c - F_{sky}
    SNR = F_{star} / \\sqrt{F_star + 2 * F_sky}

    To convert flux to magnitude, an arbitrary scale value F_0 = 1e6 e/s is used:

    M = -2.5 * log_{10}(F_{star} / F_0)

    Args:
        image (CCDData): calibrated image, pixel counts in electrons.
        stars (table-like): list of stars to be measured: auid, centroid, optional name
        r (_type_): radius of thecircular aperture
        r_in (_type_): inner radius of the annulus
        r_out (_type_): outer radius of the annulus

    Returns:
        QTable: photometry results:
                auid,
                centroid,
                flux (e/s)
                SNR (db)
                Magnitude with uncertainty,
                relative peak (max count / saturation level)

    """
    result = QTable(stars['auid', 'sky_centroid'])
    if 'name' in stars.colnames:
        result['name'] = stars['name']

    corrected = CCDData(image)

    apr = SkyCircularAperture(stars['sky_centroid'], r=r_ap)
    ap_stats = ApertureStats(corrected, apr)

    r_in, r_out = r_ann
    ann = SkyCircularAnnulus(stars['sky_centroid'],
                             r_in=r_in,
                             r_out=r_out)
    ann_stats = ApertureStats(corrected, ann)

    zero_level = 1_000_000 * image.unit / u.second
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
