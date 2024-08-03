import astropy.units as u
from astropy.table import QTable, Column, hstack
from .SimpleTransform import SimpleTransform

class BatchTransformer:
    def __init__(self, band, settings) -> None:
        self.band_ = band
        band_settings = settings[f"{self.band_[0]}{self.band_[1]}"]
        self.comp_ = band_settings['comp']
        self.check_ = band_settings['check'] if 'check' in band_settings else None

    def combine(self, provider):
        target = provider.get_target()
        batches = provider.get_batches(self.band_)
        result = [(id := batch['id'][0],
                   target[target['id'] == id][0],
                   batch[batch['auid'] == self.comp_][0],
                   None if not self.check_
                   else batch[batch['auid'] == self.check_][0],
                   SimpleTransform.create(batch, self.band_[0], self.band_[1]))
                  for batch in batches
                  ]
        return result


    def calculate(self, provider):
        data = self.combine(provider)

        transformed = [(id,
                        target['time'],
                        target['airmass'],
                        comp['auid'],
                        comp[self.band_[0]],
                        comp[self.band_[1]],
                        xfm(target, comp),
                        (xfm.Ta, xfm.Ta_err),
                        (xfm.Tab, xfm.Tab_err),
                        None if not check else check['auid'],
                        None if not check else check[self.band_[0]],
                        None if not check else check[self.band_[1]],
                        None if not check else xfm(check, comp)
                        )
                       for id, target, comp, check, xfm in data]

        def mag_err(x):
            return (x[0].value, x[1].value)

        target_result = QTable([
            Column([x[0] for x in transformed],
                   name='batch'),
            Column([x[1].jd for x in transformed],
                   name='time'),
            Column([x[2] for x in transformed],
                   name='airmass'),
            Column([x[3] for x in transformed],
                   name='comp'),
            Column([x[4] for x in transformed],
                   name=f'comp {self.band_[0]}'),
            Column([x[5] for x in transformed],
                   name=f'comp {self.band_[1]}'),
            Column([mag_err(x[6][0]) for x in transformed],
                   name=self.band_[0],
                   unit=u.mag,
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
            Column([mag_err(x[6][1]) for x in transformed],
                   name=self.band_[1],
                   unit=u.mag,
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
            Column([x[7] for x in transformed],
                   name='Ta',
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
            Column([x[8] for x in transformed],
                   name='Tab',
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
        ])

        if self.check_ is not None:
            check_result = QTable([
            Column([x[9] for x in transformed],
                   name='check'),
            Column([x[10] for x in transformed],
                   unit=u.mag,
                   name=f'check std {self.band_[0]}'),
            Column([x[11] for x in transformed],
                   unit=u.mag,
                   name=f'check std {self.band_[1]}'),
            Column([mag_err(x[12][0]) for x in transformed],
                   name=f"check {self.band_[0]}",
                   unit=u.mag,
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
            Column([mag_err(x[12][1]) for x in transformed],
                   name=f"check {self.band_[1]}",
                   unit=u.mag,
                   dtype=[('mag', 'f4'), ('err', 'f4')]),
            ])
            return hstack([target_result, check_result])
        else:
            return target_result
