import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import unittest
from unittest.mock import patch
from pathlib import Path
from vso.util import LayoutBase, ImageLayout

class LayoutBaseTest(unittest.TestCase):

    def test_constructor(self):
        root = '/home/test'

        l1 = LayoutBase(root)
        self.assertEqual(str(l1.root_), root)
        self.assertEqual(l1.create_, True)

        l2 = LayoutBase(root, create=False)
        self.assertEqual(str(l2.root_), root)
        self.assertEqual(l2.create_, False)

    def test_root_dir(self):
        root = '/home/test'

        with patch.object(Path, 'exists') as mock_exists:
            with patch.object(Path, 'mkdir') as mock_mkdir:
                mock_exists.return_value = True
                l = LayoutBase(root)
                self.assertEqual(str(l.root_dir), root)
                self.assertTrue(mock_exists.called)
                self.assertFalse(mock_mkdir.called)

        with patch.object(Path, 'exists') as mock_exists:
            with patch.object(Path, 'mkdir') as mock_mkdir:
                mock_exists.return_value = False
                l = LayoutBase(root)
                self.assertEqual(str(l.root_dir), root)
                self.assertTrue(mock_exists.called)
                self.assertTrue(mock_mkdir.called)

class ImageLayoutTest(unittest.TestCase):

    def test_constructor(self):
        root = '/home/test'
        tag='19010102'
        target='Polaris'
        l = ImageLayout(root, create=False)
        self.assertEqual(str(l.get_images(tag, target).root_dir),
                         str(Path(root) / Path(tag) / Path(target)))

if __name__ == '__main__':
    unittest.main()
