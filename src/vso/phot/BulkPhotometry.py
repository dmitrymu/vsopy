import astropy.units as u
from pathlib import Path
import ccdproc as ccdp
import concurrent.futures as cf
import itertools
import logging
from ..calibr import CalibrationMatcher
from ..reduce import load_and_solve, calibrate_image
from .measure import measure_photometry
from.BatchAggregator import BatchAggregator
from astropy.table import QTable, vstack
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

def load_settings(path):
    settings = read_json(path)
    ap = settings['aperture']
    ap_unit = u.Unit(ap['unit'])
    return Aperture(
        r=ap['r_ap'] * ap_unit,
        r_in=ap['r_in'] * ap_unit,
        r_out=ap['r_out'] * ap_unit)

def append_table(table1, table2):
    return table2 if not table1 else vstack([table1, table2])


class BulkPhotometry:
    matcher = None

    def __init__(self, session_dir, calibr_dir):
        self.session_dir_ = Path(session_dir)
        self.solved_dir_ = self.session_dir_ / 'solved'
        if not self.solved_dir_.exists():
            self.solved_dir_.mkdir(parents=True)
        self.calibr_dir_ = calibr_dir
        logging.getLogger('astropy').setLevel(logging.ERROR)
        logging.getLogger('root').setLevel(logging.ERROR)
        self.aperture_ = load_settings(self.session_dir_ / 'settings.json')
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

    def was_not_processed(self, id, file_path):
        return str(file_path) not in self.blacklist_ and id not in self.processed_

    def process_image(self, file_path, id):
        print(f"processing #{id}: {file_path}")
        image = load_and_solve(file_path, self.solved_dir_)
        cal = BulkPhotometry.matcher.match(image.header)
        reduced = calibrate_image(image,
                                  dark=cal.dark,
                                  flat=cal.flat)
        phot_table = measure_photometry(reduced, self.stars_,
                                        self.aperture_.r,
                                        self.aperture_.r_in,
                                        self.aperture_.r_out)
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
        ifc = ccdp.ImageFileCollection(image_dir)

        files = [image_dir / f for f in ifc.summary['file']]

        return [(path, executor.submit(self.process_image, path, id))
                for path, id in zip(files, counter)
                if self.was_not_processed(id, path)]


    def process(self, image_dir):
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

        aggr = BatchAggregator()
        result = aggr.aggregate(self.image_table_, self.star_table_)

        result.write(self.session_dir_ / 'photometry.ecsv')
