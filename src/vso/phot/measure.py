import astropy.units as u
import numpy as np
from .. import reduce
from ..data import CameraRegistry
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable, Column
from photutils.aperture import (ApertureStats,
                                SkyCircularAperture,
                                SkyCircularAnnulus)


def measure_photometry(image, stars, aperture, extended=False):
    """Extract aperture photometry data from the image.

    Given the calibrated image and the list of star centroids in sky coordinates,
    measure star flux using circular aperture and sky annulus.  For detailed discussion
    see Warner B.D., A Practical Guide to Lightcurve Photometry and Analysis, 4.2.

    Flux counting is performed by photutils package
    (https://photutils.readthedocs.io). Uncertainty fromthe image is propagated to
    the flux values. Flux is normalized by the image exposure to electrons per second.

    The flux for central aperture F_c contains both star and sky electrons, the flux
    for annulus F_a - sky electrons only. Number of sky electrons in the central
    aperture is

    F_{sky} = F_{sky mean} * (pixel count for central aperture)

    Then star flux in electons is

    F_{star} = F_c - F_{sky}

    and signal to noise ratio is

    SNR = F_{star} / \\sqrt{F_star + 2 * F_sky}

    To convert flux to magnitude, an arbitrary scale value F_0 = 1e6 e/s is used:

    M = -2.5 * log_{10}(F_{star} / F_0)

    Args:
        image (CCDData):            calibrated image, pixel counts in electrons.
        stars (table-like):         list of stars to be measured: AUID, centroid, optional name
        r (Quantity):               radius of the circular aperture
        r_ann (Quantity, Quantity): inner and outer radii of the annulus
        extended (bool):            return additional stats (FWHM, ellipticity etc)

    Returns:
        QTable: photometry results:
                auid,
                centroid,
                flux (e/s)
                SNR (db)
                Magnitude with uncertainty,
                relative peak (max count / saturation level)

    """
    result = QTable(stars['auid', 'radec2000'])
    if 'name' in stars.colnames:
        result['name'] = stars['name']

    corrected = CCDData(image)

    apr = SkyCircularAperture(stars['radec2000'], r=aperture.r)
    ap_stats = ApertureStats(corrected, apr)

    ann = SkyCircularAnnulus(stars['radec2000'],
                             r_in=aperture.r_in,
                             r_out=aperture.r_out)
    ann_stats = ApertureStats(corrected, ann)

    zero_level = 1_000_000 * image.unit / u.second
    exp = image.header['EXPTIME'] * u.second
    camera = CameraRegistry.get(image.header['instrume'])
    read_noise = 0 if not camera or image.unit != u.electron else camera.read_noise(image.header['GAIN'])
    Ft = ap_stats.sum
    Fb = ap_stats.sum_aper_area.value * ann_stats.mean
    FtFb = Ft - Fb
    Ft2Fb = Ft + 2 * Fb + read_noise * ap_stats.sum_aper_area.value
    result['flux'] = FtFb / exp
    # prevent NaN in snr column
    snr = np.clip(FtFb.value/np.sqrt(np.abs(Ft2Fb.value)), 1e-10, 1e30)
    result['snr'] = 10 * np.log10(snr) * u.db
    mag = -2.5 * np.log10(result['flux'] / zero_level).value

    Ft_err = ap_stats.sum_err
    Fb_err = ap_stats.sum_aper_area.value * (ann_stats.sum_err / ann_stats.sum_aper_area.value)
    flux_err = np.sqrt(Ft_err*Ft_err + Fb_err*Fb_err) / exp
    mag_err = (2.5 * flux_err / result['flux'] / np.log(10)).value
    result['M'] = Column(list(zip(mag, mag_err)),
                         unit=u.mag,
                         dtype=[('mag','f4'), ('err', 'f4')])

    result['peak'] = np.nan if not camera else ap_stats.max.value / camera.max_adu.value
    result['sky_centroid'] = ap_stats.sky_centroid

    if extended:
        result['ellipticity'] = ap_stats.ellipticity
        result['fwhm'] = ap_stats.fwhm
        result['orientation'] = ap_stats.orientation

    return result


def process_image(path, matcher, solver, centroids, aperture):
    try:
        image = reduce.update_wcs(CCDData.read(path, unit='adu'), solver(path))
        calibration = matcher.match(image.header)
        reduced = reduce.calibrate_image(image,
                                        dark=calibration.dark,
                                        flat=calibration.flat)
        return measure_photometry(reduced, centroids, aperture)
    except Exception:
        return None
