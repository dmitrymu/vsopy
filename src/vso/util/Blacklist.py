import json
from pathlib import Path

class Blacklist:
    def __init__(self, path = None):
        self.blacklist_ = {}
        if path is not None and Path(path).exists():
            with open(path) as file:
                self.blacklist_ = json.load(file)

    def contains(self, path):
        return str(path) in self.blacklist_

    def add(self, path, reason='?'):
        self.blacklist_.setdefault(str(path), str(reason))

    def save(self, path):
        with open(path, mode='w') as file:
            json.dump(self.blacklist_, file)
