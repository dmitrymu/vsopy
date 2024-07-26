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

    def fetch(self, uri):
        try:
            return open(download_file(uri, cache=self.cache_web_content))
        except Exception as ex:
            msg = f"Downloading {uri} failed: {ex}"
            return io.StringIO(msg)

    def fetch_content(self, uri):
        with self.fetch(uri) as stream:
            return stream.read()

    def _fetch_content(uri):
        def fetcher(self, *args, **kwargs):
            return self.fetch_content(uri(self, *args, **kwargs))
        return fetcher

    @_fetch_content
    def get_vsx_votable(self, star_name:str) -> str:
        return ("http://www.aavso.org/vsx/index.php?view=query.votable"
                f"&ident={name_for_api(star_name)}")

    @_fetch_content
    def get_chart_stars(self, star_name: str, fov=60, maglimit=16) -> str:
        return ("https://www.aavso.org/apps/vsp/api/chart/?format=json"
                f"&star={name_for_api(star_name)}&fov={fov}&maglimit={maglimit}")

    @_fetch_content
    def get_std_field_stars(self, ra, dec, fov=60, maglimit=16) -> str:
        return ("https://www.aavso.org/apps/vsp/api/chart/?format=json&special=std_field"
                f"&ra={ra}&dec={dec}&fov={fov}&maglimit={maglimit}")

    @_fetch_content
    def get_std_fields(self) -> str:
        return "https://www.aavso.org/vsx/index.php?view=api.std_fields&format=json"
