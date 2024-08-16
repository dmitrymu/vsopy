import json

class Settings:
    def __init__(self, path) -> None:
        self.data_ = {}
        if not path is None:
            with open(path) as file:
                self.data_ = json.load(file)

    @property
    def data(self):
        return self.data_

    @property
    def bands(self):
        return [] if 'bands' not in self.data_ else self.data_['bands']

    def get_comp(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return band_settings['comp']

    def get_check(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return  band_settings['check'] if 'check' in band_settings else None

    def set_comp(self, band, value):
        label = f"{band[0]}{band[1]}"
        # self.data_["diff_photometry"]={}
        # self.data_["diff_photometry"][label] = {}
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['comp'] = value

    def set_check(self, band, value):
        label = f"{band[0]}{band[1]}"
        # self.data_["diff_photometry"]={}
        # self.data_["diff_photometry"][label] = {}
        self.data_.setdefault("diff_photometry", {}).setdefault(label, {})['check'] = value
