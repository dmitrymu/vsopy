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

    def get_images(self, tag: str, name: str) ->Path:
        return TargetLayout(self.root_dir / Path(tag) / Path(name),
                           create=self.create_)
                           


class WorkLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def tmp_dir(self):
        return self.root_dir / 'tmp'
    
    @property
    @LayoutBase._enforce
    def calibr_dir(self):
        return self.root_dir / 'calibr'

