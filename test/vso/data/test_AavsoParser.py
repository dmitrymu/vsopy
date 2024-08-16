import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import astropy.units as u
import unittest
from vso.data import AavsoParser
from astropy.coordinates import SkyCoord, Angle

class AavsoParserTest(unittest.TestCase):

    def test_std_fields(self):
        input="""
        {"StandardFields":
            {"@version":"1",
             "StandardField":[
                {"Id":"1001","Name":"NGC 1252","RA":"47.704167","Dec":"-57.766667","Fov":"120","Count":"36"},
                {"Id":"1002","Name":"M67","RA":"132.838375","Dec":"11.800667","Fov":"20","Count":"211"},
                {"Id":"1003","Name":"NGC 3532","RA":"166.3875","Dec":"-58.73","Fov":"45","Count":"288"}
        ]}}
        """
        p = AavsoParser()
        fields = p.parse_std_fields(input)
        self.assertEqual(len(fields), 3)
        self.assertEqual(fields.colnames, ['name', 'radec2000', 'fov', 'count'])
        self.assertCountEqual(fields['name'], ['NGC 1252', 'M67', 'NGC 3532'])
        self.assertCountEqual(fields['radec2000'].ra.value,
                              [47.704167, 132.838375, 166.3875])

    def test_vsx_votable(self):
        input = """
<VOTABLE version="1.0">
<RESOURCE>
<DESCRIPTION>International Variable Star Index (VSX) Query Results</DESCRIPTION>
<TABLE>
<FIELD id="auid" name="AUID"/>
<FIELD id="name" name="Name"/>
<FIELD id="const" name="Const"/>
<FIELD id="radec2000" name="Coords(J2000)"/>
<FIELD id="varType" name="VarType"/>
<FIELD id="maxMag" name="MaxMag"/>
<FIELD id="maxPass" name="MaxMagPassband"/>
<FIELD id="minMag" name="MinMag"/>
<FIELD id="minPass" name="MinMagPassband"/>
<FIELD id="epoch" name="Epoch"/>
<FIELD id="novaYr" name="NovaYear"/>
<FIELD id="period" name="Period"/>
<FIELD id="riseDur" name="RiseDuration"/>
<FIELD id="specType" name="SpecType"/>
<FIELD id="disc" name="Discoverer"/>
<DATA>
<TABLEDATA>
<TR>
<TD>000-BDB-211</TD>
<TD>SX UMa</TD>
<TD>UMa</TD>
<TD>201.55608333,56.25697222</TD>
<TD>RRC</TD>
<TD>10.580</TD>
<TD>V</TD>
<TD>11.210</TD>
<TD>V</TD>
<TD>52746.48600</TD>
<TD/>
<TD>0.3071178</TD>
<TD>38</TD>
<TD>A4-F5</TD>
<TD>Sergei Belyavsky (1914)</TD>
</TR>
</TABLEDATA>
</DATA>
</TABLE>
</RESOURCE>
</VOTABLE>
        """
        p = AavsoParser()
        star = p.parse_vsx_votable(input)
        self.assertEqual(len(star), 1)
        self.assertEqual(star['auid'][0], '000-BDB-211')
        self.assertEqual(star['name'][0], 'SX UMa')
        self.assertEqual(star['radec2000'][0].ra.value, 201.55608333)
        self.assertEqual(star['radec2000'][0].dec.value, 56.25697222)

    def test_chart(self):
        input="""
{"chartid":"X37313LN","image_uri":"https://apps.aavso.org/vsp/chart/X37313LN.png?format=json",
 "star":"XZ Cyg","fov":60.0,"maglimit":16.0,"title":"","comment":"",
 "resolution":150,"dss":false,"special":null,
 "photometry":[
    {"auid":"000-BCG-990","ra":"19:31:18.38","dec":"56:13:11.5","label":"82",
     "bands":[{"band":"V","mag":8.204,"error":0.016},{"band":"B","mag":8.744,"error":0.035},
              {"band":"Rc","mag":7.89,"error":0.032},{"band":"Ic","mag":7.613,"error":0.038},
              {"band":"J","mag":7.14,"error":0.013},{"band":"H","mag":6.838,"error":0.015},
              {"band":"K","mag":6.817,"error":0.011}],"comments":""},
    {"auid":"000-BCH-225","ra":"19:35:47.52","dec":"55:57:33.5","label":"87",
     "bands":[{"band":"V","mag":8.685,"error":0.043},{"band":"B","mag":9.843,"error":0.061},
              {"band":"Ic","mag":7.488,"error":0.083},{"band":"J","mag":6.599,"error":0.005},
              {"band":"H","mag":6.032,"error":0.011},{"band":"K","mag":5.918,"error":0.013}],
              "comments":""},
    {"auid":"000-BMS-017","ra":"19:32:45.40","dec":"56:23:19.5","label":"109",
     "bands":[{"band":"V","mag":10.935,"error":0.029},{"band":"B","mag":11.141,"error":0.046},
              {"band":"Rc","mag":10.834,"error":0.049},{"band":"Ic","mag":10.713,"error":0.054}],
              "comments":""},
    {"auid":"000-BMS-019","ra":"19:32:29.54","dec":"56:30:14.2","label":"114",
     "bands":[{"band":"V","mag":11.401,"error":0.029},{"band":"B","mag":12.378,"error":0.042},
              {"band":"Rc","mag":10.858,"error":0.04},{"band":"Ic","mag":10.382,"error":0.045}],
              "comments":""},
    {"auid":"000-BMS-018","ra":"19:31:36.12","dec":"56:17:23.2","label":"116",
     "bands":[{"band":"V","mag":11.564,"error":0.044},{"band":"B","mag":12.092,"error":0.062},
              {"band":"Rc","mag":11.289,"error":0.08},{"band":"Ic","mag":11.011,"error":0.088}],
              "comments":""},
    {"auid":"000-BMS-020","ra":"19:33:23.18","dec":"56:18:17.1","label":"118",
     "bands":[{"band":"V","mag":11.768,"error":0.017},{"band":"B","mag":12.157,"error":0.046},
              {"band":"Rc","mag":11.532,"error":0.024},{"band":"Ic","mag":11.329,"error":0.038}],
              "comments":""}
    ],
 "auid":"000-BCH-041","ra":"19:32:29.31","dec":"56:23:17.5"}
"""
        p = AavsoParser()
        chart = p.parse_chart(input)
        self.assertEqual(len(chart), 6)
        self.assertEqual(set(chart.colnames),
                         set(['auid', 'radec2000', 'U', 'B', 'V', 'Rc', 'Ic']))
        self.assertDictEqual(chart.meta,
                             dict(chart_id='X37313LN',
                                  auid='000-BCH-041',
                                  star='XZ Cyg',
                                  radec2000=SkyCoord(
                                      ra=Angle("19:32:29.31", unit=u.hourangle),
                                      dec=Angle("56:23:17.5", unit=u.deg)
                                  )))
