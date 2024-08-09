from enum import Enum

class FrameType(Enum):
    BIAS = 'Bias'
    DARK = 'Dark'
    FLAT = 'Flat'

    def __str__(self):
        return self.value

