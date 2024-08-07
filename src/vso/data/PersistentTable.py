import numpy as np
from pathlib import Path
from astropy.table import QTable

class PersistentTable:
    """ A wrapper around QTable serialized to disk.
    """
    def __init__(self, path, initializer=QTable) ->None:
        """ Create the wrapper.

            Parameters:
            path: path-like, pointing to the serialized table representation.
            initializer: a function called with no arguments providing QTable object
                         if no serialized table exists
        """
        self.path_ = Path(path)
        self.table_ = None
        self.init_ = initializer
        self.format_ = 'ascii.ecsv'

    def flush(self) -> None:
        """ Flush the table to the disk
        """
        self.table_.write(self.path_, format=self.format_, overwrite=True)

    def get(self) -> QTable:
        """ Get underlying table object

            Result:
            table:  QTable either read from disk or provided by initializer.
        """
        if not self.table_:
            if self.path_.exists():
                self.table_ = QTable.read(self.path_, format=self.format_)
            else:
                self.table_ = self.init_()
                self.flush()
        return self.table_

    @staticmethod
    def init_from_template(template):
        """ Convenience method to create empty table with required structure.

            Parameters:
            template: QTable providing structure

            returns
            result:  Empty QTable with the same set of columns as template,
                     both names and data types
        """
        return QTable(template)[[]]

    def append(self, row):
        """ Append a row to the table and flush.

            Parameters:
            row: dictionary-like, data row to be added.

            Returns:
            table: updated QTable
        """
        self.get().add_row(row)
        self.flush()
        return self.table_

    def row_by_key(self, field, key):
        """ Get row by the value of the field.

            Parameters:
            field: the name of the table column.
            key: the value, type must match the data type of the column

            Returns:
            row:  None of table does not have table[field]==key;
                  row if key matches exactly one value in 'field'

            Raises:
                KeyError if 'key' matches multiple rows.
        """
        rows = self.table_[self.table_[field] == key]
        return None if len(rows) == 0 else rows[0] if len(rows) == 1 else KeyError(f"{len(rows)} rows found for {field}={key}")

    def row_by_keys(self, keys):
        """ Get row by the value of the field.

            Parameters:
            field: the name of the table column.
            key: the value, type must match the data type of the column

            Returns:
            row:  None of table does not have table[field]==key;
                  row if key matches exactly one value in 'field'

            Raises:
                KeyError if 'key' matches multiple rows.
        """
        filter = np.all([self.get()[field] == val for field, val in keys.items()], axis=0)
        rows = self.get()[filter]
        return None if len(rows) == 0 else rows[0] if len(rows) == 1 else KeyError(f"{len(rows)} rows found for {field}={key}")
