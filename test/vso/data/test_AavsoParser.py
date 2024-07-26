import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src')))

import unittest
from vso.data import AavsoParser

class LayoutBaseTest(unittest.TestCase):

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
