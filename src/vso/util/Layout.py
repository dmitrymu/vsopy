from pathlib import Path

class LayoutBase:
    """ Generic layout
        Supports root directory and enforces directory creation.
    """
    def __init__(self, root, create=True):
        """ Set up layout at root directory.

            root - path to root directory.
            create - whether to enforce directory creation.
        """
        self.create_ = create
        self.root_ = Path(root)

    def _enforce(dir):
        """ Decorator enforcing directory creation.

            dir - a method creating directory
        """
        def enforcer(self, *args, **kwargs):
            result = Path(dir(self, *args, **kwargs))
            if self.create_ and not result.exists():
                result.mkdir(parents=True)
            return result
        return enforcer

    @property
    @_enforce
    def root_dir(self) -> Path:
        """ Returns path to the root directory
        """
        return self.root_

class Session:
    def __init__(self, tag=None, name=None):
        self.tag_ = tag
        self.name_ = name.replace(' ', '_')

    @property
    def rel_path(self):
        return Path(self.tag_) / Path(self.name_)

    @property
    def name(self):
        return self.name_.replace('_', '')

class TargetLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def lights_dir(self):
        return self.root_dir / 'Light'

class ImageLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_images(self, session) ->Path:
        return TargetLayout(self.root_dir / session.rel_path,
                           create=self.create_)

class SessionLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def solved_dir(self):
        return self.root_dir / 'solved'

    @property
    def blacklist_file_path(self):
        return self.root_dir / 'blacklist.json'

    @property
    def chart_file_path(self):
        return self.root_dir / 'chart.ecsv'

    @property
    def centroid_file_path(self):
        return self.root_dir / 'centroids.ecsv'
    @property
    def settings_file_path(self):
        return self.root_dir / 'settings.json'

    @property
    def photometry_file_path(self):
        return self.root_dir / 'photometry.ecsv'



class WorkLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def tmp_dir(self):
        return self.root_dir / 'tmp'

    @property
    @LayoutBase._enforce
    def charts_dir(self):
        return self.root_dir / 'charts'

    @property
    @LayoutBase._enforce
    def calibr_dir(self):
        return self.root_dir / 'calibr'

    def get_session(self, session) ->Path:
        return SessionLayout(self.root_dir / Path('session') / session.rel_path,
                             create=self.create_)

