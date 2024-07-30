import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import io
import unittest
from unittest.mock import patch
from vso.data import AavsoApi

def mock_download_success(uri, **kwargs):
    return uri

def mock_open(uri, **kwargs):
    return io.StringIO(uri)

class AavsoApiTest(unittest.TestCase):

    @patch(f"vso.data.AavsoApi.download_file", mock_download_success)
    @patch(f"builtins.open", mock_open)
    def test_download(self):
        api = AavsoApi()
        self.assertEqual(
            api.get_std_fields(),
            "https://www.aavso.org/vsx/index.php?view=api.std_fields&format=json"
        )
        self.assertEqual(
            api.get_vsx_votable('SX UMa'),
            "http://www.aavso.org/vsx/index.php?view=query.votable&ident=SX+UMa"
        )
        self.assertEqual(
            api.get_star_chart('SX UMa', fov=7.5, maglimit=17.5),
            ("https://www.aavso.org/apps/vsp/api/chart/"
             "?format=json&star=SX+UMa&fov=7.5&maglimit=17.5")
        )
        self.assertEqual(
            api.get_std_field_chart(0.42, 42.0, fov=7.5, maglimit=17.5),
            ("https://www.aavso.org/apps/vsp/api/chart/"
             "?format=json&special=std_field&ra=0.42&dec=42.0&fov=7.5&maglimit=17.5")
        )
        self.assertEqual(
            api.get_chart_by_id('ABCD'),
            ("https://apps.aavso.org/vsp/api/chart/ABCD/?format=json")
        )
