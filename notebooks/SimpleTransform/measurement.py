from astropy.table import QTable, Column
from sequence import StarSequence
import numpy as np
from dataclasses import dataclass
from vsopy.util import ValErr, MagErr, MagErrDtype
from util import lognorm_params

class Measurement:
    def __init__(self, rng):
        self.rng = rng

    def __call__(self, seq: StarSequence) -> QTable:
        raise NotImplementedError


class SimpleMeasurement(Measurement):

    @dataclass
    class Params:
        offset: float = 0.0
        scatter: float = 0.01
        uncert: ValErr = ValErr(0.05, 0.05)

    def __init__(self, rng, params: Params):
        self.rng = rng
        self.params = params

    def generate_magnitudes(self, seq: StarSequence) -> QTable:
        return seq.stars['M']['mag'] + self.params.offset + \
            self.rng.normal(0, self.params.scatter, len(seq.stars))

    def generate_uncertainties(self, seq: StarSequence) -> np.ndarray:
        mu, sigma = lognorm_params(self.params.uncert.val, self.params.uncert.err)
        return self.rng.lognormal(mean=mu, sigma=sigma, size=len(seq.stars))

    def __call__(self, seq: StarSequence) -> QTable:
        measured = self.generate_magnitudes(seq)
        uncert = self.generate_uncertainties(seq)
        measured = [MagErr(m, u) for m, u in zip(measured, uncert)]
        return QTable({'auid': seq.stars['auid'],
                    'band': seq.stars['band'],
                    'instr M': Column(measured, dtype=MagErrDtype)})

class SkewedMeasurement(SimpleMeasurement):

    def __init__(self, rng, params: SimpleMeasurement.Params, skew: tuple[float, float]):
        super().__init__(rng, params)
        self.skew = skew

    def generate_magnitudes(self, seq: StarSequence) -> QTable:
        measured = super().generate_magnitudes(seq)
        ci = measured[seq.stars['band'] == seq.bands[0]] - measured[seq.stars['band'] == seq.bands[1]]
        measured[seq.stars['band'] == seq.bands[0]] += ci * self.skew[0]
        measured[seq.stars['band'] == seq.bands[1]] += ci * self.skew[1]
        return measured
