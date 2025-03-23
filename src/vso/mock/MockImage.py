import astropy.units as u
import numpy as np

from astropy import wcs
from astropy.convolution import Gaussian2DKernel
from astropy.nddata import CCDData, StdDevUncertainty
from collections import namedtuple

MockStar = namedtuple('MockStar', ['peak', 'coord', 'fwhm', 'ellipticity', 'angle'])

class MockImageBuilder:
    def __init__(self, shape, seed=42):
        self.shape_ = shape
        self.rng_ = np.random.RandomState(seed)
        self.data_ = np.zeros(self.shape_)

    @property
    def data(self):
        return self.data_

    def add_noise(self, mean, std):
        mu = np.log(mean**2 / np.sqrt(mean**2 + std**2))
        sigma = np.sqrt(np.log(1 + std**2 / mean**2))
        self.data_ += self.rng_.lognormal(mu, sigma, self.data_.size).reshape(self.shape_)

#    def add_star(self, peak, pos, fwhm, ell, angle):
    def add_star(self, star):
        sigma_x = star.fwhm / 2.355
        sigma_y = sigma_x * (1 - star.ellipticity)
        kernel =  Gaussian2DKernel( sigma_x, sigma_y, star.angle.to(u.rad))
        self.add_kernel(star.coord, star.peak, kernel)

    def add_kernel(self, pos, weight, kernel):
        w, h = self.data_.shape
        kernel_w, kernel_h = kernel.shape
        x, y = pos
        center_x, center_y = tuple(kernel.center)
        r_x = kernel_w - center_x - 1
        r_y = kernel_h - center_y - 1
        target_left = max(0, x - r_x)
        target_right = min(w-1, x + r_x)
        target_top = max(0, y - r_y)
        target_bottom = min(h-1, y + r_y)
        patch_left = max(0, r_x - x)
        patch_right = min(kernel_w-1, kernel_w - (x+1+r_x-(w-1)))
        patch_top = max(0, r_y - y)
        patch_bottom = min(kernel_h-1, kernel_h - (y+1+r_y-(h-1)))

        norm = kernel.array[kernel.center[0], kernel.center[1]]
        self.data_[target_top:target_bottom, target_left:target_right] += weight * \
            kernel.array[patch_top:patch_bottom, patch_left:patch_right] / norm

    def get_image(self, pixel_scale=1*u.arcsec, header={}, read_noise=1):
        w, h = self.shape_
        cs = wcs.WCS(naxis=2)
        cs.wcs.crpix = [w//2, h//2]
        cs.wcs.cdelt = np.array([pixel_scale.to(u.deg).value, pixel_scale.to(u.deg).value])
        cs.wcs.crval = [0, 0]
        cs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        return CCDData(self.data_,
                       wcs=cs,
                       header=header,
                       uncertainty=StdDevUncertainty(np.sqrt(self.data_ + read_noise**2),
                                                     unit=u.electron),
                       unit=u.electron)
