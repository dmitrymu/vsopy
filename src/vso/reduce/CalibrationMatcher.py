
from astropy.nddata import CCDData
from astropy.time import Time
import astropy.units as u
import ccdproc as ccdp
import numpy as np
from pathlib import Path
from collections import namedtuple
from ..util import FrameType

FRAME_LIGHT = 'Light'

MATCHER_KEYWORDS = ['instrume', 'gain', 'xbinning', 'ybinning', 'offset', 'filter',
                    'ccd-temp', 'date-obs', 'exptime', 'file']

Calibration = namedtuple('Calibration', ['bias', 'dark', 'flat'])

class FrameCollection:
    """A collection of calibration master frames ofa specific type

    Provides filtering by FITS header fields that must match
    """
    def __init__(self, calibr_dir, frame_type: str) -> None:
        """Create frame collection

        Args:
            calibr_dir ([path-like]): path to the directory containing calibration master frames.
            frame_type (str): frame type
        """
        self.dir_ = Path(calibr_dir)
        self.ifc_ = ccdp.ImageFileCollection(self.dir_).filter(frame=frame_type)

    def filter(self, header):
        """Filter master frames compatible with the header

        Args:
            header (dict-like): header of a frame to be calibrated

        Returns:
            table-like: filtered collection of files
        """
        return self.ifc_.filter(instrume = header['instrume'],
                                gain=header['gain'],
                                xbinning=header['xbinning'],
                                ybinning=header['ybinning'],
                                offset=header['offset']
                                )

class CalibrationMatcher:
    def __init__(self, calibr_dir,
                 temp_tolerance=1*u.K,
                 exposure_tolerance=.1):
        self.calibr_dir_ = Path(calibr_dir)
        self.bias_ = FrameCollection(self.calibr_dir_, FrameType.BIAS.value)
        self.dark_ = FrameCollection(self.calibr_dir_, FrameType.DARK.value)
        self.flat_ = FrameCollection(self.calibr_dir_, FrameType.FLAT.value)
        self.cache_ = {}
        self.temp_tolerance_ = temp_tolerance
        self.exposure_tolerance_ = exposure_tolerance

    def temp_filter(self, header, collection):
        dt = np.abs(collection['ccd-temp'] - header['ccd-temp']) * u.K
        dt_filter = dt <= self.temp_tolerance_
        if not np.any(dt_filter):
            raise RuntimeError(
                f"No matching frames for T = {header['ccd-temp']} C"
                f" with tolerance {self.temp_tolerance_}"
            )
        return collection[dt_filter]

    def exp_filter(self, header, collection):
        de = np.abs(collection['exptime']
                    - header['exptime']) / header['exptime']
        de_filter = de <= self.exposure_tolerance_
        if not np.any(de_filter):
            raise RuntimeError(
                f"No matching frames for exposure = {header['exptime']} sec"
                f" with tolerance {self.exposure_tolerance_:.2%}"
            )
        return collection[de_filter]

    def most_recent(self, header, collection):
        dt = Time(header['date-obs']).jd - Time(collection['date-obs']).jd
        past = collection[dt > 0]
        dt_past = dt[dt > 0]
        return past[np.argmin(dt_past)]

    def closest_time(self, header, collection):
        dt = np.abs(Time(header['date-obs']).jd - Time(collection['date-obs']).jd)
        return collection[np.argmin(dt)]

    def load_image(self, file):
        path = self.calibr_dir_ / file
        if path not in self.cache_:
            self.cache_[path] = CCDData.read(path)
        return self.cache_[path]

    def match_bias(self, header):
        candidates = self.bias_.filter(header)
        temp_filtered = self.temp_filter(header, candidates.summary)
        return self.most_recent(header, temp_filtered)['file']

    def match_dark(self, header, scale=False, future=False):
        candidates = self.dark_.filter(header)
        temp_exp_filtered = self.temp_filter(header,
                                             self.exp_filter(header, candidates.summary))
        return (self.closest_time(header, temp_exp_filtered)
                if future else
                self.most_recent(header, temp_exp_filtered))['file']

    def match_flat(self, header):
        candidates = self.flat_.filter(header).filter(filter=header['filter'])
        temp_filtered = self.temp_filter(header, candidates.summary)
        return self.most_recent(header, temp_filtered)['file']


    def match(self, header, scale=False):
        if header['frame'] == FrameType.BIAS.value:
            return Calibration(None, None, None)
        elif header['frame'] == FrameType.DARK.value:
            return Calibration(None if not scale else self.load_image(self.match_bias(header)),
                               None,
                               None)
        elif header['frame'] == FrameType.FLAT.value:
            return Calibration(None if not scale else self.load_image(self.match_bias(header)),
                               self.load_image(self.match_dark(header, scale=scale, future=True)),
                               None)
        elif header['frame'] == FRAME_LIGHT:
            return Calibration(None if not scale else self.load_image(self.match_bias(header)),
                               self.load_image(self.match_dark(header, scale=scale)),
                               self.load_image(self.match_flat(header)))
        else:
            raise RuntimeError(f"Unsupported frame type '{header['frame']}'")
