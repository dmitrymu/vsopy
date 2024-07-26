import io
from astropy.utils.data import download_file

def name_for_api(star_name):
    return star_name.replace(' ', '+')

class AavsoApi:
    """ Encapsulate HTTP APIs to AAVSO database.

        See https://www.aavso.org/apis-aavso-resources
    """

    def __init__(self, cache_web_content=True) -> None:
        self.cache_web_content = cache_web_content
        pass


    def fetch_content(self, uri):
        try:
            return open(download_file(uri, cache=self.cache_web_content))
        except Exception as ex:
            msg = f"Downloading {uri} failed: {ex}"
            return io.StringIO(msg)

    def fetch_vsx_votable(self, star_name):
        return self.fetch_content(
            "http://www.aavso.org/vsx/index.php?view=query.votable"
            f"&ident={name_for_api(star_name)}").read()

    def fetch_chart_stars(self, star_name, fov=60, maglimit=16):
        return self.fetch_content(
            "https://www.aavso.org/apps/vsp/api/chart/?format=json"
            f"&star={name_for_api(star_name)}&fov={fov}&maglimit={maglimit}").read()

    def fetch_field_stars(self, ra, dec, fov=60, maglimit=16):
        return self.fetch_content(
            "https://www.aavso.org/apps/vsp/api/chart/?format=json&special=std_field"
            f"&ra={ra}&dec={dec}&fov={fov}&maglimit={maglimit}").read()

    def fetch_std_fields(self):
        return self.fetch_content(
            "https://www.aavso.org/vsx/index.php?view=api.std_fields&format=json").read()