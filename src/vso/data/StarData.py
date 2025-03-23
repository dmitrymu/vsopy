from . import AavsoApi
from . import AavsoParser
from . import PersistentTable
import astropy.units as u
from astropy.table import QTable
from pathlib import Path


class StarData:
    """ The interface to AAVSO VSX and VSP tool

        Provides the API to download photometry data for  variable stars
        and standard fields.  Uses three levels of caching:
        - HTTP responses from AAVSO are cached by astropy.utils.data.download_file
          (see https://docs.astropy.org/en/stable/utils/data.html#cache-management)
        - Photometry tables are stored in charts directory in ECSV format
          (see https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html#astropy.io.ascii.Ecsv)
        - Photometry tables are cached in memory
    """
    def __init__(self, charts_dir, cache_web_content=True, normalize_charts=True):
        """ Create data API

            Parameters:

            charts-dir: path-like pointing to the directory to store serialized photometry data
            cache_web_content: bool, passed to astropy.utils.data.download_file
            normalize_charts: bool, return normalized chart as a tuple (centroids, sequence)
        """
        self.normalize_ = normalize_charts
        self.charts_dir_ = Path(charts_dir)
        self.api_ = AavsoApi(cache_web_content)
        self.parser_ = AavsoParser()
        self.std_fields_ = PersistentTable(
            self.charts_dir_ / 'std_fields.ecsv',
            lambda: self.parser_.parse_std_fields(
                self.api_.get_std_fields()
            )
        )
        self.charts_ = PersistentTable(
            self.charts_dir_ / 'charts.ecsv',
            lambda: PersistentTable.init_from_template(dict(
                name=[''],
                fov=[0.0]*u.arcmin,
                maglimit=[0.0]*u.mag,
                id=['']
            ))
        )
        self.charts_cache_ = {}

    @property
    def charts(self):
        """ The table containing all downloaded charts.

            Returns:
            charts: QTable containing name, FOV, magnitude limit, and AAVSO chart ID
        """
        return self.charts_.get()

    @property
    def std_fields(self):
        """ The table containing all known standard fields.

            Returns:
            std_fields: QTable containing name, center of the field coordinates,
            FOV, and star count
        """
        return self.std_fields_.get()

    def get_chart_path(self, id):
        """ Construct the path to serialized photometry data

            Returns:
            path: paths to the serialized chart in ECSV format
        """
        if self.normalize_:
            return  (self.charts_dir_ / f"{id}_c.ecsv", self.charts_dir_ / f"{id}_s.ecsv")
        else:
            return  self.charts_dir_ / f"{id}.ecsv"

    def load_chart(self, id):
        """ Access chart by ID, transparently download, cache, and serialize.

            Parameters:
            id: string, AAVSO chart ID

            Returns:
            QTable containing photometry data.
        """
        if id not in self.charts_cache_:
            if self.normalize_:
                centroids, sequence = self.parser_.parse_norm_chart(self.api_.get_chart_by_id(id))
                centr_path, seq_path = self.get_chart_path(id)
                self.charts_cache_[id] = (
                    PersistentTable(
                        centr_path,
                        lambda: centroids),
                    PersistentTable(
                        seq_path,
                        lambda: sequence)
                )
            else:
                self.charts_cache_[id] = PersistentTable(
                    self.get_chart_path(id),
                    lambda: self.parser_.parse_chart(self.api_.get_chart_by_id(id))
                )
        pt = self.charts_cache_[id]
        return (pt[0].get(), pt[1].get()) if self.normalize_ else pt.get()

    def is_std_field(self, name):
        """ Check whether the name belongs to standard field

            Returns:
            is_std_field: True if it's standard field name, False if it's variable star name.
        """
        return name in self.std_fields['name']

    def get_chart(self, name, fov=None, maglimit=16.0*u.mag):
        """ Get photometry data for either star or standard field.

            Parameters:
            name: the name of the object (star or field)
            fov: field of view in arc minutes (ignored for standard fields).
            maglimit: magnitude limit for photometry sequence

            Returns:
            chart:  QTable containing photometry data.
        """
        real_fov = self.std_fields_.row_by_key('name', name)['fov'] if self.is_std_field(name) else fov

        cached = self.charts_.row_by_keys(
                dict(name=name, fov=real_fov, maglimit=maglimit)
            )
        if not cached:
            text = None
            if self.is_std_field(name):
                field = self.std_fields_.row_by_key('name', name)
                text = self.api_.get_std_field_chart(
                    field['radec2000'].ra,
                    field['radec2000'].dec,
                    real_fov,
                    maglimit
                )
            else:
                text = self.api_.get_star_chart(name, real_fov, maglimit)

            chart = self.parser_.parse_norm_chart(text) if self.normalize_ else self.parser_.parse_chart(text)
            id = chart[1].meta['chart_id'] if self.normalize_ else chart.meta['chart_id']
            self.charts_.append(dict(
                name=name,
                fov=real_fov,
                maglimit=maglimit,
                id=id
            ))
            if self.normalize_:
                centr_path, seq_path = self.get_chart_path(id)
                self.charts_cache_[id] = (
                    PersistentTable(
                        centr_path,
                        lambda: chart[0]),
                    PersistentTable(
                        seq_path,
                        lambda: chart[1])
                )
            else:
                self.charts_cache_[id] = PersistentTable(
                    self.get_chart_path(id),
                    lambda: chart
                )
            pt = self.charts_cache_[id]
            return (pt[0].get(), pt[1].get()) if self.normalize_ else pt.get()
        else:
            return self.load_chart(cached['id'])


