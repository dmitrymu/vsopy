import unittest
import astropy.units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord
from unittest.mock import patch
from vso.data import StarData, PersistentTable, AavsoApi, AavsoParser
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src')))



STD_FIELDS = QTable(dict(
    name=['SA00'],
    radec2000=SkyCoord(ra=[0]*u.deg, dec=[0]*u.deg),
    fov=[60.0]*u.arcmin,
    count=[0]
))


CHARTS_TEMPLATE = dict(
    name=[''],
    fov=[0.0],
    maglimit=[0.0],
    id=['']
)

STAR_CHART = """
{"chartid":"X37313LN",
 "photometry":[
    {"auid":"000-BMS-020","ra":"19:33:23.18","dec":"56:18:17.1","label":"118",
     "bands":[{"band":"V","mag":11.768,"error":0.017},{"band":"B","mag":12.157,"error":0.046},
              {"band":"Rc","mag":11.532,"error":0.024},{"band":"Ic","mag":11.329,"error":0.038}],
              "comments":""}
 ]
}
"""

class StarDataTest(unittest.TestCase):

    @patch.object(PersistentTable, "get")
    def test_construct_from_initializer(self, mock_get):
        mock_get.side_effect = [PersistentTable.init_from_template(CHARTS_TEMPLATE),
                                STD_FIELDS, STD_FIELDS
                                ]
        sd = StarData('/home/test')
        charts = sd.charts
        mock_get.assert_called_once()
        self.assertEqual(len(charts), 0)
        self.assertEqual(set(charts.colnames), set(['name', 'fov', 'maglimit', 'id']))
        fields = sd.std_fields
        self.assertEqual(len(fields), 1)
        self.assertEqual(set(fields.colnames), set(['name', 'fov', 'radec2000', 'count']))
        self.assertTrue(sd.is_std_field('SA00'))

    @patch(f"vso.data.PersistentTable.QTable.write")
    @patch.object(AavsoApi, 'get_star_chart')
    @patch.object(AavsoParser, 'parse_std_fields')
    def test_get_star_chart(self, parse_std_fields, mock_get_star_chart, mock_write):
        parse_std_fields.return_value = STD_FIELDS
        mock_get_star_chart.return_value = STAR_CHART
        sd = StarData('/home/test')
        chart = sd.get_chart('A Star', 60*u.arcmin, 15*u.mag)
        self.assertEqual(len(chart), 1)

    @patch(f"vso.data.PersistentTable.QTable.write")
    @patch.object(AavsoApi, 'get_std_field_chart')
    @patch.object(AavsoParser, 'parse_std_fields')
    def test_get_star_chart(self, parse_std_fields, mock_get_std_field_chart, mock_write):
        parse_std_fields.return_value = STD_FIELDS
        mock_get_std_field_chart.return_value = STAR_CHART
        sd = StarData('/home/test', normalize_charts=False)
        chart = sd.get_chart('SA00')
        self.assertEqual(len(chart), 1)

    @patch(f"vso.data.PersistentTable.QTable.write")
    @patch.object(AavsoApi, 'get_std_field_chart')
    @patch.object(AavsoParser, 'parse_std_fields')
    def test_get_norm_star_chart(self, parse_std_fields, mock_get_std_field_chart, mock_write):
        parse_std_fields.return_value = STD_FIELDS
        mock_get_std_field_chart.return_value = STAR_CHART
        sd = StarData('/home/test')
        centroids, sequence = sd.get_chart('SA00')
        self.assertEqual(len(centroids), 1)
        self.assertEqual(len(sequence), 4)
