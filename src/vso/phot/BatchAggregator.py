import numpy as np
import numpy.lib.recfunctions as rf
from astropy.table import join, unique, join_distance
from astropy.time import Time, TimeDelta
from ..util import ordered_bands

def column_mask(column):
    return column.mask if hasattr(column, 'mask') else np.full((len(column)), False)

class BatchAggregator:
    """Prepare per-image photometry for differential photometry .

    The basic assumption is that filer is changed for each new image
    (e.g. B V R B V R B V R), so images are naturally grouped into multi-band
    batches (e,g. BVR, BVR, BVR) and every batch is a single observation
    for differential photometry. The order of filter change must be consistent
    inside a whole session.
    """
    def __init__(self) -> None:
        """Create the processor

        """

    def batch(self, bands, image_table):
        """Combine single-band images into multi-band batches.

        Args:
            bands (iterable of str): list of bands of interest
            image_table (_type_): _description_

        Returns:
            _type_: _description_
        """
        image_table['finish'] = image_table['time'].jd + TimeDelta(image_table['exposure']).jd
        image_table.add_index('band')

        per_band = {}
        exposure = {}
        for band in bands:
            images = image_table.loc[band]
            per_band[band] = images[f'id', f'time', f'finish', f'airmass']
            per_band[band]['time'] = per_band[band]['time'].jd
            exposure[band] = np.max(images['exposure'].value)
            to_rename = [f'id', f'finish', f'airmass']
            per_band[band].rename_columns(to_rename, [f"{c}_{band}" for c in to_rename])

        batches = None
        prev_band = None

        for band in bands:
            batches = per_band[band] if not batches else join(
                batches, per_band[band],
                keys='time',
                table_names=[prev_band, band],
                join_funcs={'time': join_distance(
                    (exposure[prev_band] + exposure[band])/86400
                )},
                join_type='outer'
            )
            if prev_band is not None:
                batches.remove_column('time_id')
                mask = (column_mask(batches[f"time_{band}"])
                        | column_mask(batches[f"time_{prev_band}"]))
                batches = batches[~mask]
                batches.rename_column(f"time_{band}", 'time')
            prev_band = band

        batches.rename_column('time', f"time_{bands[-1]}")
        start_time = np.min(rf.structured_to_unstructured(batches[[f"time_{b}" for b in bands]].as_array()), axis=1)
        finish_time = np.max(rf.structured_to_unstructured(batches[[f"finish_{b}" for b in bands]].as_array()), axis=1)
        batches['time'] = Time((start_time + finish_time) / 2, format='jd')
        batches['duration'] = TimeDelta(finish_time - start_time).to('second')
        batches['airmass'] = np.mean(rf.structured_to_unstructured(batches[[f"airmass_{b}" for b in bands]].as_array()), axis=1)
        batches.remove_columns([f"time_{b}" for b in bands])
        batches.remove_columns([f"finish_{b}" for b in bands])
        batches['id'] = range(1, len(batches) + 1)
        return batches

    def aggregate(self, image_table, star_table):
        """Create batches and join per-star photometry


        Args:
            image_table (table-like): per-image info
                                      (band, exposure, gain, observation time, air mass)
            star_table (table-like): per-star photometry
                                     (auid, band, magnitude, SNR, peak)

        Returns:
            QTable: batched photometry
                    (auid, band, mag B, SNR B, peak B, mag V, SNR V, ...)
        """
        bands = ordered_bands(image_table['band'])
        batches = self.batch(bands, image_table)


        stargroup =  join(batches['id','time', 'duration', 'airmass'],
                          unique(star_table['auid',]),
                          join_type='cartesian')
        selector = stargroup
        for band in bands:
            selector = join(selector,
                    batches['id', f"id_{band}"],
                    keys='id')

        result = selector
        for band in bands:
            s = star_table[star_table['band'] == band]['id', 'auid', 'M', 'snr', 'peak']
            to_rename = ['id', 'snr', 'peak']
            s.rename_columns(to_rename, [f"{c}_{band}" for c in to_rename])
            s.rename_column('M', f"instr {band}")
            result = join(result, s, keys=[f"id_{band}", 'auid'], join_type='inner')

        result.remove_columns([f"id_{band}" for band in bands])
        result.rename_columns([c for c in result.colnames if '_' in c], [c.replace('_', ' ') for c in result.colnames if '_' in c])
        return result

