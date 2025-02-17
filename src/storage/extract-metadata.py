import logging
import sys
import sqlite3

from astropy.io import fits
from astropy.time import Time
from pathlib import Path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger('default')

STORAGE_ROOT=Path('/mnt/astro/source')
FILE_ROOT=Path('/mnt/public/astro/archive')
DATABASE = Path('/srv/public/images')

class DictTable:
    def __init__(self, db, table):
        self.db_ = db
        self.cur_ = db.cursor()
        self.table_ = table
        self.data_ = dict({name: id for id, name in self.cur_.execute(
            f'SELECT id, name FROM {self.table_};'
        )})

    def id(self, name):
        if name not in self.data_:
            self.cur_.execute(
                f'INSERT INTO {self.table_}(name) VALUES (?)',
                (name,)
            )
            self.data_[name] = self.cur_.lastrowid
            self.db_.commit()
        return self.data_[name]


class Processor:

    def __init__(self, db):
        self.db_ = sqlite3.connect(db)
        self.cur_ = self.db_.cursor()
        self.cameras_ = DictTable(self.db_, 'camera')
        self.objects_ = DictTable(self.db_, 'object')
        self.frame_types_ = DictTable(self.db_, 'frame_type')


    def folders(self, parent):
        if not parent:
            return self.cur_.execute(
                'SELECT id, path FROM storage_folder WHERE parent IS NULL;'
            )
        else:
            return self.cur_.execute(
                'SELECT id, path FROM storage_folder WHERE parent = ?;',
                (parent,)
            )

    def insert_dir(self, parent, dir):
        self.cur_.execute(
            'INSERT INTO storage_folder(parent, path) VALUES (?, ?)',
            (parent, str(dir))
        )
        id = self.cur_.lastrowid
        logger.info(f'Inserted {id}: {dir}')
        self.db_.commit()
        return id

    def image_attr(self, path, dir, file):
        with fits.open(path) as hdul:
            header = hdul[0].header
            return (
                dir,
                str(file),
                float(Time(header['DATE-OBS']).jd),
                self.objects_.id( header.get('OBJECT', 'None')),
                self.cameras_.id(header['INSTRUME']),
                self.frame_types_.id(header['FRAME']),
                header.get('FILTER', None),
                header['EXPTIME'],
                header['GAIN'],
                header['OFFSET'],
                header['CCD-TEMP'],
                header.get('AIRMASS', None)
            )

    def insert_image(self, attr):
        self.cur_.execute(
            'INSERT INTO image(folder, file, time, object, camera, frame_type'
            '                  filter, exposure, gain, offset, temperature, airmass)'
            ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            attr
        )

    def update_image(self, id, attr):
        self.cur_.execute(
            'UPDATE image SET folder = ?, file = ?, time = ?, object = ?, camera = ?,'
            '                 frame_type = ?, filter = ?, exposure = ?, gain = ?,'
            '                 offset = ?, temperature = ?, airmass = ?'
            ' WHERE id = ?',
            attr + (id,)
        )

    def process_image(self, dir, path):
        file = str(path.name)
        logger.info(f'processing {path}')
        attributes = self.image_attr(FILE_ROOT / path.relative_to(STORAGE_ROOT), dir, file)
        self.cur_.execute(f"SELECT id FROM image WHERE file = ? AND folder = ?", (file, dir))
        id = self.cur_.fetchone()[0]
        if id is None:
            self.insert_image(attributes)
        else:
            self.update_image(id, attributes)
        self.db_.commit()

    def traverse(self, root, parent=None):
        logger.info(f'processing {root}')
        dirs = {path: id for (id, path) in self.folders(parent)}
        for d in root.iterdir():
            if d.is_dir():
                name = str(d.relative_to(STORAGE_ROOT))
                id = self.insert_dir(parent, name) if name not in dirs else dirs[name]
                self.traverse(d, id)
            else:
                self.process_image(parent, d)


def main():
    processor = Processor(DATABASE)
    processor.traverse(STORAGE_ROOT)

if __name__ == '__main__':
    sys.exit(main())
