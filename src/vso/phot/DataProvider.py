import astropy.units as u
import numpy as np
from astropy.table import join, vstack, QTable, Column
from ..util import ordered_bands

def no_nan_mag(batch, band_a, band_b):
    filter_nan = np.all([~np.isnan(batch[band_a]['mag']),
                         ~np.isnan(batch[band_a]['err']),
                         ~np.isnan(batch[band_b]['mag']),
                         ~np.isnan(batch[band_b]['err'])], axis=0)
    return batch[filter_nan]

class Filter:
    def __init__(self, comparator) -> None:
        self.comparator_ = comparator

    def __call__(self, data, band):
        return np.all([self.comparator_(data, b) for b in band], axis=0)

    @staticmethod
    def saturation(threshold):
        return Filter(lambda data, b: data[f"peak {b}"] < threshold)

    @staticmethod
    def snr(threshold):
        return Filter(lambda data, b: data[f"snr {b}"] > threshold)

    @staticmethod
    def time(min, max):
        return Filter(lambda data, b: (data["time"] > min) & (data["time"] < max))

    @staticmethod
    def passthrough():
        return Filter(lambda data, *args: np.full([len(data)], True))

class DataProvider:
    """ Prepare photometry data for transformation.
    """

    def __init__(self, photometry, chart,
                 bands = None,
                 filters = [
                     Filter.saturation(.9),
                     Filter.snr(10.0 * u.db)
                 ],
                 tag = None):
        self.bands_ = ordered_bands(bands)
        self.photometry_ = photometry
        self.chart_ = chart
        self.data_ = None if photometry is None or chart is None else join(
            self.photometry_, chart[self.bands_ + ['auid']], keys='auid')
        if self.data_ is not None and tag is not None:
            self.data_["tag"] = [tag] * len(self.data_)
        self.filters_ = filters if len(filters) > 0 else [Filter.passthrough()]

    def get_batches(self, bands):
        data = self.data_ if self.filters_ is None else self.data_[
            np.logical_and.reduce([f(self.data_, bands) for f in self.filters_])
        ]
        return no_nan_mag(data, bands[0], bands[1]).group_by('id').groups

    @property
    def target_name(self):
        return self.chart_.meta['star']

    @property
    def chart_id(self):
        return self.chart_.meta['chart_id']

    @property
    def bands(self):
        return self.bands_

    def get_target(self):
        return self.photometry_[self.photometry_['auid'] == self.chart_.meta['auid']]

    @staticmethod
    def combine(coll):
        def offset_id(data, offset):
            result = QTable(data)
            result['id'] += offset
            return result

        result = DataProvider(None, None)
        bands = set()
        filters = []
        for provider in coll:
            result.data_ = provider.data_ if result.data_ is None else vstack([
                result.data_, offset_id(provider.data_, np.max(result.data_['id']))
            ])
            bands.update(provider.bands_)
            filters.extend(provider.filters_)
        result.bands_ = ordered_bands(bands)
        result.filters_ = filters
        return result

