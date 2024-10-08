import astropy.units as u
from pathlib import Path
import ccdproc as ccdp
import concurrent.futures as cf
import itertools
import logging
from ..reduce import load_and_solve, calibrate_image, CalibrationMatcher
from ..util import ordered_bands
from .measure import measure_photometry
from.BatchAggregator import BatchAggregator
from astropy.table import QTable, vstack, unique
import psutil
import json
from collections import namedtuple
from astropy.time import Time

Aperture = namedtuple('Aperture', ['r', 'r_in', 'r_out'])


def num_workers():
    GIGABYTE = 1024*1024*1024
    vm = psutil.virtual_memory()
    return 2 * vm.available // (GIGABYTE*3)

def read_json(path):
    with open(path) as file:
       return json.load(file)

def read_table_if_exists(path):
    return QTable.read(path) if path.exists() else None

def append_table(table1, table2):
    return table2 if not table1 else vstack([table1, table2])


class BulkPhotometry:
    """Extract photometry data from images in two-level image directory.

    Image directory is assumed to contain per-band subdirectories.
    Calibration images are stored in a separate directory. Parameters and
    results are placed in a separate directory.

    Source image is calibrated and plate solved before flux measurement.
    Command line version of ASTAP (https://www.hnsky.org/astap.htm) is
    used as plate solver. Plate solving results are cached in a
    subdirectory of session directory to speed up repeated processing.

    Image processing is parallelized, currently one process per 1.5 GB
    of free memory. (TODO: pass number of process to constructor).

    Results of the process are three tables stored in the session directory.
    * stars.ecsv - per-star photometry data
    * images - per-image data like air mass and observation time
    * photometry - per-band photometry grouped into multi-band observations
                   (see BatchAggregator for details).

    The actual output is the third file that is the input for differential
    photometry algorithms.

    First two files are intermediate results. If processing is run repeatedly,
    already processed files are skipped and new results are appended.

    Files that fail processing are marked in black list stored in the session
    directory and skipped if processing repeats.
    """
    matcher = None

    def __init__(self, session_dir, calibr_dir):
        self.session_dir_ = Path(session_dir)
        self.solved_dir_ = self.session_dir_ / 'solved'
        if not self.solved_dir_.exists():
            self.solved_dir_.mkdir(parents=True)
        self.calibr_dir_ = calibr_dir
        logging.getLogger('astropy').setLevel(logging.ERROR)
        logging.getLogger('root').setLevel(logging.ERROR)
        self.aperture_ = self.load_aperture()
        self.stars_ = QTable.read(self.session_dir_ / 'centroids.ecsv')
        self.star_table_ = read_table_if_exists(self.star_table_path)
        self.image_table_ = read_table_if_exists(self.image_table_path)
        self.blacklist_ = read_json(self.blacklist_path) if self.blacklist_path.exists() else {}
        self.processed_ = set() if not self.image_table_ else set(self.image_table_['id'])

    @property
    def star_table_path(self):
        return self.session_dir_ / 'stars.ecsv'

    @property
    def image_table_path(self):
        return self.session_dir_ / 'images.ecsv'

    @property
    def blacklist_path(self):
        return self.session_dir_ / 'blacklist.json'

    @staticmethod
    def worker_init(me):
        BulkPhotometry.matcher = CalibrationMatcher(me.calibr_dir_,
                                                    temp_tolerance=2*u.K)

    def load_aperture(self):
        settings = read_json(self.session_dir_ / 'settings.json')
        ap = settings['aperture']
        ap_unit = u.Unit(ap['unit'])
        return Aperture(
            r=ap['r_ap'] * ap_unit,
            r_in=ap['r_in'] * ap_unit,
            r_out=ap['r_out'] * ap_unit)

    def update_bands(self, bands):
        path = self.session_dir_ / 'settings.json'
        settings = read_json(path)
        settings['bands'] = bands
        with open(path, mode='w') as f:
            json.dump(settings, f, indent=2)

    def was_not_processed(self, id, file_path):
        return str(file_path) not in self.blacklist_ and id not in self.processed_

    def process_image(self, file_path, id: int):
        """Process a single image.

        Load, plate solve, calibrate, and measure photometry.

        Args:
            file_path (path-like): path to the image
            id (int): image ID.

        Returns:
            tuple: First element is QTable containing per-star photometry data
                   Second element is a dictionary  containing the image properties
                   (id, path, band, exposure, gain, observation time, air mass).
        """
        print(f"processing #{id}: {file_path}")
        image = load_and_solve(file_path, self.solved_dir_)
        cal = BulkPhotometry.matcher.match(image.header)
        reduced = calibrate_image(image,
                                  dark=cal.dark,
                                  flat=cal.flat)
        phot_table = measure_photometry(reduced, self.stars_,
                                        self.aperture_.r,
                                        (self.aperture_.r_in, self.aperture_.r_out))
        band = image.header['filter']
        phot_table['id'] = [id]*len(phot_table)
        phot_table['band'] = [band]*len(phot_table)
        file_info = {
            'id': id,
            'path': str(file_path),
            'band': band,
            'exposure': image.header['exptime'] * u.second,
            'gain': image.header['gain'],
            'time': Time(image.header['date-obs']),
            'airmass': image.header['airmass']
        }
        return (phot_table, file_info)

    def process_directory(self, executor, image_dir, counter):
        """Process images from the directory

        Args:
            executor (concurrent.futures.Executor): the interface to submit
                                                    a concurrent task .
            image_dir (path-like): the directory being processed
            counter (iterable): sequence of image IDs (ID is unique
                              per session)

        Returns:
            iterable containing pairs (path to file,
                                       Future to get processing results)
        """
        ifc = ccdp.ImageFileCollection(image_dir)

        files = [image_dir / f for f in ifc.summary['file']]

        return [(path, executor.submit(self.process_image, path, id))
                for path, id in zip(files, counter)
                if self.was_not_processed(id, path)]


    def process(self, image_dir) -> None:
        """Process images from subdirectories of image directory.

        The result is stored in session directory. Images are processed
        in parallel using ProcessPoolExecutor (https://docs.python.org/3/library/concurrent.futures.html)
        If the process crashes and kills its parent terminal, the reason
        is lack of free RAM.

        Args:
            image_dir (path-like): _description_
        """
        nw = num_workers()
        counter = itertools.count(1)
        with cf.ProcessPoolExecutor(initializer=BulkPhotometry.worker_init,
                                    initargs=(self,),
                                    max_workers=nw) as executor:
            subdirs = [self.process_directory(executor,
                                              d,
                                              counter)
                       for d in image_dir.iterdir()
                       if d.is_dir()]

        def get_result(file_result):
            path = None
            try:
                path, future = file_result
                return future.result()
            except Exception as e:
                self.blacklist_[str(path)] = str(e)
                return None

        result = [x for x in [get_result(r)
                              for r in itertools.chain(*subdirs)]
                  if x is not None]

        with open(self.blacklist_path, mode="w+") as file:
            json.dump(self.blacklist_, file, indent=2)

        zipped_result = tuple(zip(*result))
        if  len(zipped_result) > 0:
            phot_list, info_list = zipped_result

            self.star_table_ = append_table(self.star_table_, vstack(phot_list))
            self.image_table_ = append_table(self.image_table_, QTable(info_list))
            self.star_table_.write(self.star_table_path, format='ascii.ecsv', overwrite=True)
            self.image_table_.write(self.image_table_path, format='ascii.ecsv', overwrite=True)

        self.update_bands(ordered_bands(self.star_table_['band']) if len(self.star_table_) > 0 else [])

        aggr = BatchAggregator()
        result = aggr.aggregate(self.image_table_, self.star_table_)

        result.write(self.session_dir_ / 'photometry.ecsv')
