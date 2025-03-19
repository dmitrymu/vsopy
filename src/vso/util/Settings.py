import astropy.units as u
import json

def convert_to_unit(x, unit):
    return x.to(unit) if hasattr(x, 'unit') else x * unit

class Aperture:
    def __init__(self, r, r_in, r_out, unit=u.arcsec):
        self.r_ = convert_to_unit(r, unit)
        self.r_in_ = convert_to_unit(r_in, unit)
        self.r_out_ = convert_to_unit(r_out, unit)

    @property
    def r(self):
        return self.r_

    @property
    def r_in(self):
        return self.r_in_

    @property
    def r_out(self):
        return self.r_out_

    @staticmethod
    def fromdict(d):
        return Aperture(d['r_ap'], d['r_in'], d['r_out'], u.Unit(d['unit']))

    def todict(self):
        unit = self.r_.unit
        return dict(r_ap=float(self.r.to(unit).value),
                    r_in=float(self.r_in.to(unit).value),
                    r_out=float(self.r_out.to(unit).value),
                    unit=str(unit))

class Settings:
    def __init__(self, path) -> None:
        self.data_ = {}
        if not path is None:
            with open(path) as file:
                self.data_ = json.load(file)
        self.disabled_ = set(self.data_['disabled']) if 'disabled' in self.data_ else set()

    @property
    def data(self):
        return self.data_

    @property
    def aperture(self):
        ap = self.data_['aperture']
        return Aperture.fromdict(ap)

    @property
    def bands(self):
        return [] if 'bands' not in self.data_ else self.data_['bands']

    def get_comp(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return band_settings['comp']

    def get_check(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return  band_settings['check'] if 'check' in band_settings else None

    def set_aperture(self, aperture):
        self.data_.setdefault('aperture', {}).update(aperture.todict())

    def set_comp(self, band, value):
        label = f"{band[0]}{band[1]}"
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['comp'] = value

    def set_check(self, band, value):
        label = f"{band[0]}{band[1]}"
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['check'] = value

    def is_star_enabled(self, star):
        return star not in self.disabled_

    def disable_star(self, star):
        self.disabled_.add(star)
