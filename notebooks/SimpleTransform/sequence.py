from astropy.table import QTable, Column, vstack
from vsopy.phot import create_simple_transform
from vsopy.util import MagErr, MagErrDtype, ValErr
from util import lognorm_params
import numpy as np
from dataclasses import dataclass


class StarSequence:

    def __init__(self,
                 stars,
                 bands
                 ):
        self.stars = stars
        self.bands = bands
        seq_A = self.stars[self.stars['band'] == self.bands[0]]
        self.comp_auid = seq_A[np.argmin(seq_A['M']['mag'])]['auid']
        self.comp_star = self.stars[self.stars['auid'] == self.comp_auid]
        self.ensemble = self.stars[self.stars['auid'] != self.comp_auid]

