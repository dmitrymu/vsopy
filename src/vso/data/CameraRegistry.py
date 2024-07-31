import astropy.units as u
import numpy as np

def ccd_gain(score, full_well, adc_bits):
    scale = (1 << adc_bits) * u.adu
    return full_well / scale / (10 ** (score / (20 * u.db)))

UNITY_XFM = lambda g, full_well, adc_bits: 1.0*u.electron/u.adu

ASI_XFM = lambda g, full_well, adc_bits : ccd_gain(g * u.db / 10, full_well, adc_bits)

class Camera:
    def __init__(self, full_well, adc_bits, xfm=UNITY_XFM, read_noise=None):
        self.full_well_ = full_well
        self.adc_bits_ = adc_bits
        self.xfm_ = xfm
        self.read_noise_ = read_noise

    def gain_to_e(self, score):
        return self.xfm_(score, self.full_well_, self.adc_bits_)

    def read_noise(self, score):
        return np.interp(score, self.read_noise_[0], self.read_noise_[1]) \
               if self.read_noise_ else None

    @property
    def adu_scale(self):
        return 1 << (16 - self.adc_bits_)

    @property
    def max_adu(self):
        return (1 << self.adc_bits_) * u.adu

class CameraRegistry:
    KNOWN_CAMERAS = {
        'ZWO CCD ASI533MM Pro': Camera(50_000 * u.electron,
                                       14,
                                       ASI_XFM,
                                       read_noise=([0, 50, 99, 100, 200, 400],
                                                   [3.8, 3.4, 3.0, 1.46, 1.27, 0.9] * u.electron))
    }
    def __init__(self):
        pass

    @staticmethod
    def get(name):
        return CameraRegistry.KNOWN_CAMERAS[name] \
               if name in CameraRegistry.KNOWN_CAMERAS \
               else None
