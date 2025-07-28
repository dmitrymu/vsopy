import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import unittest

import astropy.units as u
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import QTable
from vsopy.mock import MockImageBuilder, MockStar
from vsopy.phot import measure_photometry
from vsopy.util import Aperture

class MeasurePhotometryTest(unittest.TestCase):

    def test_one_star(self):
        SHAPE = (31,31)
        NOISE_MEAN = 100
        NOISE_STDDEV = 10
        STAR_PEAK = 1000
        STAR_POS = (15, 15)
        STAR_FWHM = 3
        STAR_AUID = 'mock-star-1'

        builder = MockImageBuilder(SHAPE)
        builder.add_noise(NOISE_MEAN, NOISE_STDDEV)
        builder.add_star(MockStar(STAR_PEAK, STAR_POS, STAR_FWHM, 0, 0*u.deg))
        image = builder.get_image(1.2 * u.arcsec,
                        dict(
                            EXPTIME=2,
                            GAIN=100,
                            instrume="ZWO CCD ASI533MM Pro"
                        ))
        centroids = QTable(dict(
            auid = [STAR_AUID],
            radec2000 = SkyCoord(ra=[0] * u.arcsec, dec=[-0] * u.arcsec)
        ))

        ph = measure_photometry(image, centroids, Aperture(5, 10, 15))

        self.assertSequenceEqual(
            ph.colnames, ['auid', 'radec2000', 'flux', 'snr', 'M', 'peak', 'sky_centroid'])
        self.assertEqual(len(ph), 1)
        self.assertEqual(ph['auid'][0], STAR_AUID)
