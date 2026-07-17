from astropy.table import QTable, Column, vstack
from sequence import StarSequence
import numpy as np
from dataclasses import dataclass
from vsopy.util import ValErr, MagErr, MagErrDtype
from vsopy.phot import SimpleTransform, create_simple_transform
from util import lognorm_params
from measurement import Measurement

class Experiment:

    @dataclass
    class GenParams:
        num_stars: int
        mag: tuple[float, float]
        ci: tuple[float, float]
        bands: tuple[str, str] = ('A', 'B')
        uncert: ValErr = ValErr(0.05, 0.05)

    def __init__(self, seed):
        self.seed = seed
        self.rng = np.random.default_rng(seed)

    def generate_sequence(self, params: GenParams) -> StarSequence:
        magnitudes = self.rng.uniform(params.mag[0], params.mag[1], params.num_stars)
        color_indices = self.rng.uniform(params.ci[0], params.ci[1], params.num_stars)
        mu, sigma = lognorm_params(params.uncert.val, params.uncert.err)
        uncert_A = self.rng.lognormal(mean=mu, sigma=sigma, size=params.num_stars)
        uncert_B = self.rng.lognormal(mean=mu, sigma=sigma, size=params.num_stars)
        auid = [f"star-{str(i).zfill(3)}" for i in range(params.num_stars)]

        def mag_column(mags, uncerts):
            return [MagErr(m, u) for m, u in zip(mags, uncerts)]

        A = QTable({'auid': auid,
                    'band': Column([params.bands[0]] * params.num_stars),
                    'M': Column(mag_column(magnitudes, uncert_A), dtype=MagErrDtype)})
        B = QTable({'auid': auid,
                    'band': Column([params.bands[1]] * params.num_stars),
                    'M': Column(mag_column(magnitudes + color_indices, uncert_B), dtype=MagErrDtype)})
        result = vstack([A, B])
        result.sort(['auid', 'band'])
        return StarSequence(result, params.bands)

    def create_target(self,
                      mag:tuple[float, float],
                      uncert:tuple[float, float],
                      bands: tuple[str, str]) -> QTable:
        a = QTable({'auid': ['target'],
                    'band': Column([bands[0]] ),
                    'M': Column([MagErr(mag[0], uncert[0])], dtype=MagErrDtype)})
        b = QTable({'auid': ['target'],
                    'band': Column([bands[1]] ),
                    'M': Column([MagErr(mag[1], uncert[1])], dtype=MagErrDtype)})
        result = vstack([a, b])
        result.sort(['auid', 'band'])
        return StarSequence(result, bands)

    def create_measurement(self, class_obj, **kwargs) -> Measurement:
        return class_obj(self.rng, **kwargs)

    def measure_ensemble(self, seq: StarSequence, meas: Measurement) -> tuple[QTable, QTable]:
        # return (ensemble, comp)
        measured = meas(seq)
        return (measured[measured['auid'] != seq.comp_auid],
                measured[measured['auid'] == seq.comp_auid])

    def create_transforms(self,
                          seq: StarSequence,
                          meas: Measurement,
                          num_measurements: int) -> list[SimpleTransform]:
        measured_ensembles = [self.measure_ensemble(seq, meas)[0] for _ in range(num_measurements)]
        bands = seq.bands
        std_A = seq.ensemble[seq.ensemble['band'] == bands[0]]['M']['mag']
        std_B = seq.ensemble[seq.ensemble['band'] == bands[1]]['M']['mag']
        return [create_simple_transform(std_A, std_B,
                                s[s['band'] == bands[0]]['instr M']['mag'],
                                s[s['band'] == bands[1]]['instr M']['mag']) for s in measured_ensembles]
