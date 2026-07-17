from dataclasses import dataclass
import itertools
from util import simple_transform_stats
from astropy.table import QTable
from experiment import Experiment, Measurement

class SequenceParamGrid(Experiment):

    @dataclass
    class GridParams:
        num_stars: list[int]
        mag_range: list[float]
        ci_range: list[float]
        base_mag: float = 10.0
        base_ci: float = 0.0

    @dataclass
    class ParamSet:
        num_stars: int
        mag: tuple[float, float]
        ci: tuple[float, float]

    def __init__(self, seed, grid_params: GridParams,
                #  meas_params: Experiment.MeasurementParams,
                 num_attempts: int,
                 num_measurements:int):
        super().__init__(seed)
        # self.meas_params = meas_params
        self.num_attempts = num_attempts
        self.num_measurements = num_measurements
        self.grid = [SequenceParamGrid.ParamSet(
            n,
            (grid_params.base_mag,
             grid_params.base_mag + mag_range),
            (grid_params.base_ci,
             grid_params.base_ci + ci_range))
            for n in grid_params.num_stars
            for mag_range in grid_params.mag_range
            for ci_range in grid_params.ci_range]


    def run_cell(self, params, meas: Measurement) -> dict:
        xfm = list(itertools.chain.from_iterable(
            [self.create_transforms(
                self.generate_sequence(
                Experiment.GenParams(
                params.num_stars,
                params.mag,
                params.ci)),
                meas,
                self.num_measurements) for _ in range(self.num_attempts)]))
        return (simple_transform_stats(xfm) |
                {'Nstars': params.num_stars,
                'Mag_range': params.mag[1] - params.mag[0],
                'CI_range': params.ci[1] - params.ci[0]})

    def run(self, meas: Measurement) -> QTable:
        return QTable([self.run_cell(params, meas) for params in self.grid])