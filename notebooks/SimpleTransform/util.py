import numpy as np
from vsopy.phot import SimpleTransform

def lognorm_params(mean: float, stddev: float) -> tuple[float, float]:
    mu = np.log((mean ** 2) / np.sqrt(stddev ** 2 + mean ** 2))
    sigma = np.sqrt(np.log(stddev ** 2 / (mean ** 2) + 1))
    return mu, sigma

def simple_transform_stats(xfm: list[SimpleTransform]) -> dict[str, float]:
    return {'Tab': np.mean([x.Tab.val for x in xfm]),
            'Tab.std': np.std([x.Tab.val for x in xfm]),
            'Tab.err': np.mean([x.Tab.err for x in xfm]),
            'Ta': np.mean([x.Ta.val for x in xfm]),
            'Ta.std': np.std([x.Ta.val for x in xfm]),
            'Ta.err': np.mean([x.Ta.err for x in xfm]),
            'Tb': np.mean([x.Tb.val for x in xfm]),
            'Tb.std': np.std([x.Tb.val for x in xfm]),
            'Tb.err': np.mean([x.Tb.err for x in xfm])}

