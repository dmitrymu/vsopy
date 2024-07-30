from . import AavsoApi
from . import AavsoParser
from . import PersistentTable
import astropy.units as u
from astropy.table import QTable
from pathlib import Path


class StarData:
    def __init__(self, charts_dir, cache_web_content=True):
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
        return self.charts_.get()

    @property
    def std_fields(self):
        return self.std_fields_.get()

    def get_chart_path(self, id):
        return  self.charts_dir_ / f"{id}.ecsv"

    def load_chart(self, id):
        if id not in self.charts_cache_:
            self.charts_cache_[id] = PersistentTable(
                self.get_chart_path(id),
                lambda: self.parser_.parse_chart(
                    self.api_.get_chart_by_id(id)
                )
            )
        return self.charts_cache_[id].get()

    def is_std_field(self, name):
        return name in self.std_fields['name']

    def get_std_field_dict(self, name, maglimit):
        field = self.std_fields_.row_by_key('name', name)
        print( field)
        return dict(ra=field['radec2000'].ra,
                    dec=field['radec2000'].dec,
                    fov=field['fov'],
                    maglimit=maglimit)

    def get_chart(self, name, fov=None, maglimit=15.0*u.mag):
        if name not in self.charts['name']:
            text = None
            real_fov = fov
            if self.is_std_field(name):
                field = self.std_fields_.row_by_key('name', name)
                real_fov = field['fov']
                text = self.api_.get_std_field_chart(
                    field['radec2000'].ra,
                    field['radec2000'].dec,
                    real_fov,
                    maglimit
                )
            else:
                text = self.api_.get_star_chart(name, real_fov, maglimit)

            chart = self.parser_.parse_chart(text)
            id = chart.meta['chart_id']
            self.charts_.append(dict(
                name=name,
                fov=real_fov,
                maglimit=maglimit,
                id=id
            ))
            self.charts_cache_[id] = PersistentTable(
                self.get_chart_path(id),
                lambda: chart
            )
            return self.charts_cache_[id].get()
        else:
            return self.load_chart(self.charts_.row_by_key('name', name)['id'])


