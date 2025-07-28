import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import unittest
from unittest.mock import patch
from pathlib import Path
from vsopy.util import LayoutBase, ImageLayout, Session, WorkLayout, SessionLayout

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

class SessionTest(unittest.TestCase):

    def test_constructor(self):
        tag='19010102'
        target='Polaris'
        session = Session(tag=tag, name=target)

        self.assertEqual(str(session.rel_path),
                         str(Path(tag) / Path(target)))

class ImageLayoutTest(unittest.TestCase):

    def test_layout(self):
        root = '/home/test'
        tag='19010102'
        target='Polaris'
        session = Session(tag=tag, name=target)

        l = ImageLayout(root, create=False)
        self.assertEqual(str(l.get_images(session).root_dir),
                         str(Path(root) / Path(tag) / Path(target)))

class WorkLayoutTest(unittest.TestCase):

    def test_layout(self):
        root = '/home/test'
        tag='19010102'
        target='Polaris'
        session = Session(tag=tag, name=target)
        l = WorkLayout(root, create=False)

        self.assertEqual(str(l.root_dir), str(Path(root)))
        self.assertEqual(str(l.calibr_dir), str(Path(root) / 'calibr'))
        self.assertEqual(str(l.charts_dir), str(Path(root) / 'charts'))

        self.assertEqual(str(l.get_session(session).root_dir),
                         str(Path(root) / 'session' / Path(tag) / Path(target)))

class SessionLayoutTest(unittest.TestCase):

    def test_layout(self):
        root = '/home/test/session'
        tag='19010102'
        target='Polaris'
        session = Session(tag=tag, name=target)
        l = SessionLayout(Path(root) / session.rel_path, create=False)

        self.assertEqual(str(l.root_dir),
                         str(Path(root) / Path(tag) / Path(target)))
        self.assertEqual(str(l.solved_dir), str(l.root_dir / 'solved'))
        self.assertEqual(str(l.blacklist_file_path), str(l.root_dir / 'blacklist.json'))
        self.assertEqual(str(l.batches_file_path), str(l.root_dir / 'batches.ecsv'))
        self.assertEqual(str(l.batch_images_file_path), str(l.root_dir / 'batch_images.ecsv'))
        self.assertEqual(str(l.settings_file_path), str(l.root_dir / 'settings.json'))
        self.assertEqual(str(l.centroid_file_path), str(l.root_dir / 'centroids.ecsv'))
        self.assertEqual(str(l.chart_file_path), str(l.root_dir / 'chart.ecsv'))
        self.assertEqual(str(l.sequence_file_path), str(l.root_dir / 'sequence.ecsv'))
        self.assertEqual(str(l.images_file_path), str(l.root_dir / 'images.ecsv'))
        self.assertEqual(str(l.measured_file_path), str(l.root_dir / 'measured.ecsv'))
        self.assertEqual(str(l.photometry_file_path), str(l.root_dir / 'photometry.ecsv'))

if __name__ == '__main__':
    unittest.main()
