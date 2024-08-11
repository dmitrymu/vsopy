
DEFAULT_BAND_ORDER = ['U', 'B', 'V', 'Rc', 'Ic']

def ordered_bands(bands):
    return (DEFAULT_BAND_ORDER
            if bands is None else
            [x for x in DEFAULT_BAND_ORDER if x in set(bands)])

def band_pairs(bands):
    return list(zip(bands[0:-1], bands[1:]))
