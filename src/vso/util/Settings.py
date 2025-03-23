import astropy.units as u
import json

def convert_to_unit(x, unit, equiv=None):
    return (x.to(unit) if equiv is None else x.to(unit, equiv))  if hasattr(x, 'unit') else x * unit

class Aperture:
    def __init__(self, r, r_in, r_out, unit=u.arcsec, scale=None):
        self.r_ = convert_to_unit(r, unit, scale)
        self.r_in_ = convert_to_unit(r_in, unit, scale)
        self.r_out_ = convert_to_unit(r_out, unit, scale)

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
    def from_dict(d):
        return Aperture(d['r_ap'], d['r_in'], d['r_out'], u.Unit(d['unit']))

    def to_dict(self):
        unit = self.r_.unit
        return dict(r_ap=float(self.r.to(unit).value),
                    r_in=float(self.r_in.to(unit).value),
                    r_out=float(self.r_out.to(unit).value),
                    unit=str(unit))

    def to_pixels(self, scale):
        return Aperture(self.r_, self.r_in_, self.r_out_,
                        u.pixel, u.pixel_scale(scale))

class PhotometrySettings:
    def __init__(self, data):
        self.data_ = data

    @property
    def comp(self):
        return self.data_['comp']

    @property
    def check(self):
        return  self.data_['check'] if 'check' in self.data_ else None

    @property
    def start(self):
        return self.data_['start']

    @property
    def finish(self):
        return self.data_['finish']

    def set_comp(self, value):
        self.data_.setdefault("comp", value)

    def set_check(self, value):
        self.data_.setdefault("check", value)

    def set_start(self, value):
        self.data_["start"] = value

    def set_finish(self, value):
        self.data_["finish"] = value

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
        return Aperture.from_dict(ap)

    @property
    def bands(self):
        return [] if 'bands' not in self.data_ else self.data_['bands']

    def photometry(self, band):
        label = f"{band[0]}{band[1]}"
        return PhotometrySettings(self.data_.setdefault("diff_photometry", {}).setdefault(label, {}))

    def get_comp(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return band_settings['comp']

    def get_check(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return  band_settings['check'] if 'check' in band_settings else None

    def set_aperture(self, aperture):
        self.data_.setdefault('aperture', {}).update(aperture.to_dict())

    def set_comp(self, band, value):
        label = f"{band[0]}{band[1]}"
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['comp'] = value

    def set_check(self, band, value):
        label = f"{band[0]}{band[1]}"
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['check'] = value

    def is_star_enabled(self, star):
        return (star) not in self.disabled_

    def disable_star(self, star):
        self.disabled_.add(star)
        self.data_['disabled'] = list(self.disabled_)

    def enable_star(self, star):
        if star in self.disabled_:
            self.disabled_.remove(star)
            self.data_['disabled'] = list(self.disabled_)
