import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src')))

import unittest
from unittest.mock import patch, Mock
from vso.calibr import CalibrationMatcher
from astropy.table import Table
from collections import namedtuple

MOCK_BIAS = Table({
    'file':['fb1'],
    'frame':['Bias'],
    'instrume':['Cam1'],
    'gain':[100],
    'xbinning':[1],
    'ybinning':[1],
    'offset':[10],
    'filter':[''],
    'ccd-temp': [-9.9],
    'date-obs': ['2024-07-20T07:00:00'],
    'exptime': [.000032]
})

MOCK_DARK = Table({
    'file':['fd1'],
    'frame':['Dark'],
    'instrume':['Cam1'],
    'gain':[100],
    'xbinning':[1],
    'ybinning':[1],
    'offset':[10],
    'filter':[''],
    'ccd-temp': [-10.0],
    'date-obs': ['2024-07-20T07:05:00'],
    'exptime': [10]
})

MOCK_FLAT = Table({
    'file':['ff1'],
    'frame':['Flat'],
    'instrume':['Cam1'],
    'gain':[100],
    'xbinning':[1],
    'ybinning':[1],
    'offset':[10],
    'filter':['V'],
    'ccd-temp': [-10.1],
    'date-obs': ['2024-07-20T07:10:00'],
    'exptime': [1.5]
})

MockIfc = namedtuple('MockIfc', ['summary'])

class CalibrationMatcherTest(unittest.TestCase):

    @patch('vso.calibr.CalibrationMatcher.CCDData.read')
    @patch('vso.calibr.CalibrationMatcher.FrameCollection')
    def test_construct(self, mock_frames, mock_read):

        colls=[Mock(summary=MOCK_BIAS), Mock(summary=MOCK_DARK), Mock(summary=MOCK_FLAT)]
        colls[0].filter.return_value=colls[0]
        colls[1].filter.return_value=colls[1]
        colls[2].filter.return_value=colls[2]
        frames = [Mock(), Mock(), Mock()]
        frames[0].filter.return_value = colls[0]
        frames[1].filter.return_value = colls[1]
        frames[2].filter.return_value = colls[2]
        mock_frames.side_effect = frames
        mock_read.return_value = []
        m = CalibrationMatcher('home/test')
        c = m.match({
            'frame':'Light',
            'instrume':'Cam1',
            'gain':100,
            'xbinning':1,
            'ybinning':1,
            'offset':10,
            'filter':'',
            'ccd-temp': -10.2,
            'date-obs': '2024-07-20T08:00:00',
            'exptime': 10
        })
        self.assertIsNone(c.bias)
        self.assertIsNotNone(c.dark)
        self.assertIsNotNone(c.flat)
