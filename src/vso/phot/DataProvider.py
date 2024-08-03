import astropy.units as u
import numpy as np
from astropy.table import join

def no_nan_mag(batch, band_a, band_b):
    filter_nan = np.all([~np.isnan(batch[band_a]['mag']),
                         ~np.isnan(batch[band_a]['err']),
                         ~np.isnan(batch[band_b]['mag']),
                         ~np.isnan(batch[band_b]['err'])], axis=0)
    return batch[filter_nan]

class DataProvider:
    """ Prepare photometry data for transformation.
    """

    def __init__(self, photometry, chart,
                 saturation_threshold=.9,
                 snr_threshold=10.0 * u.db,
                 bands=['B', 'V', 'Rc', 'Ic']):
        self.bands_ = bands
        self.photometry_ = photometry
        self.chart_ = chart
        self.data = join(self.photometry_, chart[self.bands_ + ['auid']], keys='auid')
        self.saturation_threshold = saturation_threshold
        self.snr_threshold=snr_threshold

    def get_batches(self, bands):
        saturation_filter = np.all([self.data[f"peak {b}"] < self.saturation_threshold
                                    for b in bands], axis=0)
        snr_filter = np.all([self.data[f"snr {b}"] > self.snr_threshold
                             for b in bands], axis=0)
        return no_nan_mag(self.data[saturation_filter & snr_filter], bands[0], bands[1]).group_by('id').groups

    @property
    def target_name(self):
        return self.chart_.meta['star']

    @property
    def chart_id(self):
        return self.chart_.meta['chart_id']

    def get_target(self):
        return self.photometry_[self.photometry_['auid'] == self.chart_.meta['auid']]
