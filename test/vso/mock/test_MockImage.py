import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import unittest

import astropy.units as u
import numpy as np
from vso.mock import MockImageBuilder, MockStar

class MockImageBuilderTest(unittest.TestCase):

    def test_construction(self):
        SHAPE = (31,31)
        builder = MockImageBuilder(SHAPE)
        image = builder.get_image()
        self.assertIsNotNone(image)
        self.assertEqual(len(image.header), 0)
        self.assertEqual(image.data.shape, SHAPE)
        self.assertTrue(np.all(image.data == 0))
        self.assertTrue(np.all(image.uncertainty.array == 1))
        self.assertEqual(image.wcs.naxis, 2)
        np.testing.assert_equal(image.wcs.wcs.crval, [0., 0.])
        np.testing.assert_equal(image.wcs.wcs.crpix, [15., 15.])

    def test_add_noise(self):
        SHAPE = (31,31)
        NOISE_MEAN = 100
        NOISE_STDDEV = 10
        builder = MockImageBuilder(SHAPE)
        builder.add_noise(NOISE_MEAN, NOISE_STDDEV)
        image = builder.get_image()
        np.testing.assert_almost_equal(np.mean(image.data), NOISE_MEAN, decimal=0)
        np.testing.assert_almost_equal(np.std(image.data), NOISE_STDDEV, decimal=1)
        np.testing.assert_almost_equal(np.mean(image.uncertainty.array), np.sqrt(NOISE_MEAN), decimal=0)

    def test_add_star(self):
        SHAPE = (31,31)
        STAR_PEAK = 1000
        STAR_POS = (15, 15)
        STAR_FWHM = 3
        builder = MockImageBuilder(SHAPE)
        builder.add_star(MockStar(STAR_PEAK, STAR_POS, STAR_FWHM, 0, 0*u.deg))
        image = builder.get_image()
        np.testing.assert_almost_equal(image.data[STAR_POS[0], STAR_POS[1]], STAR_PEAK)
        np.testing.assert_almost_equal(image.data[0, 0], 0)
