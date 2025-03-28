import astropy.units as u
import json
import numpy as np
import xml.etree.ElementTree as ET
from astropy.coordinates import SkyCoord, Angle
from astropy.table import QTable, Column
from typing import Tuple

VSX_VOTABLE_FIELDS = set([
    'auid', 'name', 'const', 'radec2000', 'varType',
    'maxMag', 'maxPass', 'minMag', 'minPass'
])

def extract_metadata(chart):
    meta=dict(
        chart_id=chart['chartid'],
    )
    if 'auid' in chart:
        meta['auid'] = chart['auid']
    if 'star' in chart:
        meta['star'] = chart['star']
    if 'ra' in chart and 'dec' in chart:
        meta['radec2000'] = SkyCoord(
            ra=Angle(chart['ra'], unit=u.hourangle),
            dec=Angle(chart['dec'], unit=u.deg)
        )
    return meta



class AavsoParser:
    """ Parse data received from AAVSO HTTP APIs
    """

    def __init__(self) -> None:
        pass

    def parse_std_fields(self, text: str) -> QTable:
        """Parse standard fields list

        Args:
            text (str): JSON text returned by AavsoApi.get_std_fields

        Returns:
            QTable: List of standard fields
        """
        data = json.loads(text)

        fields = [f for f in data['StandardFields']['StandardField']]
        result = dict(
            name=[f['Name'] for f in fields],
            radec2000=SkyCoord(ra=[Angle(f['RA'], unit=u.deg) for f in fields],
                               dec=[Angle(f['Dec'], unit=u.deg) for f in fields]),
            fov=[float(f['Fov']) * u.arcmin for f in fields],
            count=[int(f['Count']) for f in fields]
        )
        return QTable(result)

    def parse_vsx_votable(self, xml: str) -> QTable:
        """Parse star data from AAVSO VSX

        Args:
            xml (str): XML VOTABLE text from AavsoApi.get_vsx_votable

        Raises:
            RuntimeError: if VOTABLE is empty or contains multiple strings

        Returns:
            QTable: Single-row table containing AAVSO AUID, star name and coordinates,
                and other data
        """
        root = ET.fromstring(xml)
        table = root.find('RESOURCE').find('TABLE')
        rows = list(table.find('DATA').find('TABLEDATA').iter('TR'))
        num_rows = len(rows)
        if num_rows != 1:
            raise RuntimeError(
                f"Expected one row in VSX VOTABLE, found {num_rows}")

        fields = [(name.attrib['id'], value.text)
                  for name, value in zip(table.iter('FIELD'),
                                         rows[0].iter('TD'))]

        def parse_coord(text):
            tokens = text.split(',')
            return SkyCoord(Angle(tokens[0], unit=u.deg),
                            Angle(tokens[1], unit=u.deg))

        result = {name:
                  [parse_coord(value)] if name == 'radec2000'
                   else [float(value)]*u.mag if name in ['maxMag', 'minMag']
                      else [value]
                  for name, value in fields
                  if name in VSX_VOTABLE_FIELDS}
        return QTable(result)

    def parse_chart(self, text: str,
                    band_set=set(['U', 'B', 'V', 'Rc', 'Ic'])) -> QTable:
        """Parse chart photometry data from AAVSO VSP

        Args:
            text (str): JSON text returned by AavsoApi.get_star_chart,
                        AavsoApi.get_std_field_chart, or AavsoApi.get_chart_by_id
            band_set ([str], optional): list of bands to be stored.
                   Defaults to set(['U', 'B', 'V', 'Rc', 'Ic']).

        Raises:
            RuntimeError: if received data indicates any error.

        Returns:
            QTable: photometry data
        """
        chart = json.loads(text)
        if 'errors' in chart:
            raise RuntimeError(f"Error from AAVSO API for target {meta['name']}: "
                            f"{';'.join(chart['errors'])}")
        if len(chart['photometry']) == 0:
            return None

        [(star['auid'],
          Angle(star['ra'], unit=u.hourangle),
          Angle(star['dec'], unit=u.deg),
          )
         for star in chart['photometry']]

        band_values = [{band['band']: (float(band['mag']), float(band['error']))
                        for band in star['bands']}
                       for star in chart['photometry']]

        bands = {band:
                 Column([star[band] if band in star else (np.nan, np.nan) for star in band_values],
                        unit=u.mag,
                        dtype=[('mag', 'f4'), ('err', 'f4')])
                 for band in band_set}

        # TODO: explore masked columns

        result = dict(
            auid=[star['auid'] for star in chart['photometry']],
            radec2000=SkyCoord(
                ra=[Angle(star['ra'], unit=u.hourangle)
                    for star in chart['photometry']],
                dec=[Angle(star['dec'], unit=u.deg)
                     for star in chart['photometry']]
            )
        )

        result.update(bands)

        return QTable(result, meta=extract_metadata(chart))

    def parse_norm_chart(self, text: str) -> Tuple[QTable, QTable]:
        """Parse chart photometry data from AAVSO VSP to normalized form

        Args:
            text (str): JSON text returned by AavsoApi.get_star_chart,
                        AavsoApi.get_std_field_chart, or AavsoApi.get_chart_by_id

        Raises:
            RuntimeError: if received data indicates any error.

        Returns:
            Tuple of two QTable:
            - 'centroids': auid, RA, Dec
            - 'sequence': auid, band, magnitude, error
        """
        chart = json.loads(text)
        if 'errors' in chart:
            raise RuntimeError(f"Error from AAVSO API for target {meta['name']}: "
                            f"{';'.join(chart['errors'])}")
        if len(chart['photometry']) == 0:
            return (None, None)

        centroids = dict(
            auid=[star['auid'] for star in chart['photometry']],
            radec2000=SkyCoord(
                ra=[Angle(star['ra'], unit=u.hourangle)
                    for star in chart['photometry']],
                dec=[Angle(star['dec'], unit=u.deg)
                     for star in chart['photometry']]
            )
        )

        bands = [(star, band) for star in chart['photometry'] for band in star['bands']]
        sequence = dict(
            auid = [star['auid'] for star, _ in bands],
            band = [band['band'] for _, band in bands],
            M = Column([(float(band['mag']), float(band['error'])) for _, band in bands],
                        unit=u.mag,
                        dtype=[('mag', 'f4'), ('err', 'f4')])
        )

        return QTable(centroids), QTable(sequence, meta=extract_metadata(chart))
