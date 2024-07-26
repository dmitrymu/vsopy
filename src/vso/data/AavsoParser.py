import astropy.units as u
import json
import xml.etree.ElementTree as ET
from astropy.coordinates import SkyCoord, Angle
from astropy.table import QTable

class AavsoParser:
    """ Parse data received from AAVSO HTTP APIs
    """
    VSX_VOTABLE_FIELDS = set([
        'auid', 'name', 'const', 'radec2000', 'varType',
        'maxMag', 'maxPass', 'minMag', 'minPass'
        ])


    def __init__(self) -> None:
        pass

    def parse_std_fields(self, text):
        """ Parse standard fields list.

            Expects JSON text returned by AavsoApi.fetch_std_fields
        """
        data = json.loads(text)

        fields = [f for f in data['StandardFields']['StandardField']]
        result = dict(
            name = [f['Name'] for f in fields],
            radec2000 = SkyCoord(ra=[Angle(f['RA'], unit=u.deg) for f in fields],
                                 dec=[Angle(f['Dec'], unit=u.deg) for f in fields]),
            fov  = [float(f['Fov']) * u.arcmin for f in fields],
            count = [int(f['Count']) for f in fields]
            )
        return QTable(result)
    

    def parse_vsx_votable(self, xml):
        """ Parse star data from AAVSO VSX

            Expects XML VOTABLE text.
        """
        root = ET.fromstring(xml)
        table = root.find('RESOURCE').find('TABLE')
        rows = list(table.find('DATA').find('TABLEDATA').iter('TR'))
        num_rows = len(rows)
        if num_rows != 1:
            raise Exception(f"Expected one row in VSX VOTABLE, found {num_rows}")
        
        fields = [(name.attrib['id'], value.text)
            
                  for name, value in zip(table.iter('FIELD'),
                                         rows[0].iter('TD'))]
        
        def parse_coord(text):
            tokens = text.split(',')
            return SkyCoord(Angle(tokens[0], unit=u.deg),
                            Angle(tokens[1], unit=u.deg))
        
        result = {name : 
                  [parse_coord(value) if name == 'radec2000'
                    else float(value) if name in ['maxMag', 'minMag']
                      else value]
                  for name, value in fields
                  if name in AavsoParser.VSX_VOTABLE_FIELDS}
        return QTable(result)
