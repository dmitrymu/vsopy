import astropy.units as u
import ccdproc as ccdp
import numpy as np
from .. import util
from ..data import CameraRegistry
from astropy.nddata import CCDData


def calibrate_image(image, bias=None, dark=None, flat=None):

    camera_name = image.header['instrume']
    camera = CameraRegistry.get(camera_name)
    adu_scale = 1 if camera is None else camera.adu_scale

    source = CCDData(image)
    source.data = source.data / adu_scale

    image_gain = image.header['gain']
    e_gain = camera.gain_to_e(image_gain) if camera else None
    e_noise = camera.read_noise(image_gain) if camera else None

    reduced = ccdp.ccd_process(source,
                               error=True,
                               gain=e_gain,
                               readnoise=e_noise,
                               exposure_key='exptime',
                               exposure_unit=u.second,
                               dark_frame=dark,
                               master_bias=bias,
                               master_flat=flat)

    reduced.meta['bias-sub'] = 'T' if bias else 'F'
    reduced.meta['dark-sub'] = 'T' if dark else 'F'
    reduced.meta['flat-sub'] = 'T' if flat else 'F'
    reduced.meta['comment'] = 'Created by VSO reduction pipeline'
    return reduced
