from vso.data import PersistentTable
from unittest.mock import patch
from astropy.table import QTable
import numpy as np
import unittest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src')))

TABLE_PATH = 'home/test/table.ecsv'

def mock_qtable_read(path, **kwargs):
    raise NotImplementedError('QTable.read')

def mock_qtable_write(self, path, **kwargs):
    raise NotImplementedError('QTable.write')

class PersistentTableTest(unittest.TestCase):

    @patch(f"vso.data.PersistentTable.Path.exists")
    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_construct_from_initializer(self, mock_write, mock_exists):
        mock_exists.return_value = False
        t = PersistentTable(TABLE_PATH)
        self.assertEqual(str(t.path_), TABLE_PATH)
        table = t.get()
        self.assertEqual(len(table), 0)
        self.assertEqual(len(table.columns), 0)
        mock_exists.assert_called_once()
        mock_write.assert_called_once_with(t.path_, format='ascii.ecsv', overwrite=True)

    @patch(f"vso.data.PersistentTable.Path.exists")
    @patch(f"vso.data.PersistentTable.QTable.read")
    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_construct_from_persistent(self, mock_write, mock_read, mock_exists):
        mock_exists.return_value = True
        t = PersistentTable(TABLE_PATH)
        self.assertEqual(str(t.path_), TABLE_PATH)
        mock_read.return_value = QTable(dict(
            col_a = [42]
        ))
        table = t.get()
        self.assertEqual(len(table), 1)
        self.assertEqual(len(table.columns), 1)
        self.assertEqual(table['col_a'], [42])
        mock_exists.assert_called_once()
        mock_read.assert_called_once_with(t.path_, format='ascii.ecsv')
        mock_write.assert_not_called()

    @patch(f"vso.data.PersistentTable.Path.exists")
    @patch(f"vso.data.PersistentTable.QTable.read")
    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_get(self, mock_write, mock_read, mock_exists):
        mock_exists.return_value = True
        t = PersistentTable(TABLE_PATH)
        self.assertEqual(str(t.path_), TABLE_PATH)
        mock_read.return_value = QTable(dict(
            col_a = [42]
        ))
        table = t.get()
        self.assertEqual(len(table), 1)
        self.assertEqual(len(table.columns), 1)
        self.assertEqual(table['col_a'], [42])
        mock_exists.assert_called_once()
        mock_read.assert_called_once_with(t.path_, format='ascii.ecsv')
        mock_write.assert_not_called()
        table = t.get()
        self.assertEqual(len(table), 1)
        self.assertEqual(len(table.columns), 1)
        self.assertEqual(table['col_a'], [42])

    def test_init_from_template(self):
        template = dict(
            a=[''],
            b=[0],
            c=[0.0]
        )
        table = PersistentTable.init_from_template(template)
        self.assertEqual(len(table), 0)
        self.assertEqual(set(template.keys()), set(table.colnames))
        self.assertEqual(table['a'].dtype, np.dtype('<U1'))
        self.assertEqual(table['b'].dtype, np.dtype(int))
        self.assertEqual(table['c'].dtype, np.dtype(float))

    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_append_empty(self, mock_write):
        row = dict(
            a=np.array(['']),
            b=np.array([0]),
            c=np.array([0.0])
        )
        table = PersistentTable(TABLE_PATH,
                                initializer=lambda: PersistentTable.init_from_template(row))
        # t = table.get()
        # mock_write.assert_called_once()
        # self.assertEqual(len(t), 0)
        # self.assertEqual(set(row.keys()), set(t.colnames))
        t = table.append(row)
        self.assertEqual(mock_write.call_count, 2)
        self.assertEqual(len(t), 1)
        self.assertEqual(dict(t[0]), row)


    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_append(self, mock_write):
        row = dict(
            a=np.array(['']),
            b=np.array([0]),
            c=np.array([0.0])
        )
        table = PersistentTable(TABLE_PATH,
                                initializer=lambda: PersistentTable.init_from_template(row))
        t = table.get()
        mock_write.assert_called_once()
        self.assertEqual(len(t), 0)
        self.assertEqual(set(row.keys()), set(t.colnames))
        t = table.append(row)
        self.assertEqual(mock_write.call_count, 3)
        self.assertEqual(len(t), 1)
        self.assertEqual(dict(t[0]), row)

    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_row_by_key(self, _):
        template = dict(
            name=[''],
            id=[0],
            value=[0]
        )
        row1 = dict(
            name='a',
            id=1,
            value=42
        )
        row2 = dict(
            name='b',
            id=2,
            value=4242
        )
        table = PersistentTable(TABLE_PATH,
                                initializer=lambda: PersistentTable.init_from_template(template))
        t = table.get()
        self.assertEqual(len(t), 0)
        t = table.append(row1)
        t = table.append(row2)
        self.assertEqual(len(t), 2)
        self.assertDictEqual(dict(table.row_by_key('name', 'a')), row1)
        self.assertDictEqual(dict(table.row_by_key('id', 2)), row2)
        self.assertIsNone(table.row_by_key('id', 42))

    @patch(f"vso.data.PersistentTable.QTable.write")
    def test_row_by_keys(self, _):
        template = dict(
            name=[''],
            id=[0],
            value=[0]
        )
        row1 = dict(
            name='a',
            id=1,
            value=42
        )
        row2 = dict(
            name='a',
            id=3,
            value=44
        )
        row3 = dict(
            name='b',
            id=2,
            value=4242
        )
        table = PersistentTable(TABLE_PATH,
                                initializer=lambda: PersistentTable.init_from_template(template))
        t = table.get()
        self.assertEqual(len(t), 0)
        t = table.append(row1)
        t = table.append(row2)
        t = table.append(row3)
        self.assertEqual(len(t), 3)
        self.assertDictEqual(dict(table.row_by_keys(dict(name='a', value=42))), row1)
        self.assertIsNone(table.row_by_keys(dict(name='a', value=40)))
