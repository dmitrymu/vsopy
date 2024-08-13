import json

class Settings:
    def __init__(self, path) -> None:
        self.data_ = {}
        if not path is None:
            with open(path) as file:
                self.data_ = json.load(file)

    @property
    def bands(self):
        return [] if 'bands' not in self.data_ else self.data_['bands']

    def get_comp(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return band_settings['comp']

    def get_check(self, band):
        band_settings = self.data_["diff_photometry"][f"{band[0]}{band[1]}"]
        return  band_settings['check'] if 'check' in band_settings else None
