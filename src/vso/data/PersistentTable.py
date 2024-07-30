from pathlib import Path
from astropy.table import QTable

class PersistentTable:
    def __init__(self, path, initializer=QTable):
        self.path_ = Path(path)
        self.table_ = None
        self.init_ = initializer
        self.format_ = 'ascii.ecsv'

    def flush(self):
        self.table_.write(self.path_, format=self.format_, overwrite=True)

    def get(self) -> QTable:
        if not self.table_:
            if self.path_.exists():
                self.table_ = QTable.read(self.path_, format=self.format_)
            else:
                self.table_ = self.init_()
                self.flush()
        return self.table_

    @staticmethod
    def init_from_template(template):
        return QTable(template)[[]]

    def append(self, row):
        self.get().add_row(row)
        self.flush()
        return self.table_

    def row_by_key(self, field, key):
        rows = self.table_[self.table_[field] == key]
        return None if len(rows) == 0 else rows[0] if len(rows) == 1 else Exception(f"{len(rows)} rows found for {field}={key}")
