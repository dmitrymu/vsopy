from astropy.table import QTable, join

class BatchDataProvider:
    def __init__(self, session):
        self.batches_ = QTable.read(session.batches_file_path)
        self.batch_images_ = QTable.read(session.batch_images_file_path)
        self.images_ = QTable.read(session.images_file_path)
        self.measured_ = QTable.read(session.measured_file_path)
        self.sequence_ = QTable.read(session.sequence_file_path)

    @property
    def target_auid(self):
        return self.sequence_.meta['auid'] if 'auid' in self.sequence_.meta else None

    def instr(self, band):
        return f"instr {band}"

    def targ(self, band):
        return f"targ {band}"

    def sequence_band(self, band):
        M = self.sequence_[self.sequence_['band'] == band]['auid', 'M']
        M.rename_column('M', band)
        return M

    def sequence_band_pair(self, band_pair):
        bandA, bandB = band_pair
        return join(self.sequence_band(bandA), self.sequence_band(bandB), 'auid')

    def check_band_pair(self, band_pair, auid):
        bandA, bandB = band_pair
        s = join(self.sequence_band(bandA), self.sequence_band(bandB), 'auid')
        return s[s['auid'] == auid][bandA, bandB]

    def batch_band(self, band):
        fltr = self.images_[self.images_['filter'] == band]['image_id',]
        image_M = join(fltr, self.measured_['image_id', 'auid', 'M'])
        batch_M = join(image_M, self.batch_images_, 'image_id')['batch_id', 'auid', 'M']
        batch_M.rename_column('M', self.instr(band))
        return batch_M

    def batch_band_pair(self, band_pair):
        bandA, bandB = band_pair
        return join(self.batch_band(bandA), self.batch_band(bandB), ['batch_id', 'auid'])

    def batch_and_sequence_band_pair(self, band_pair):
        return join(self.batch_band_pair(band_pair), self.sequence_band_pair(band_pair), 'auid')

    def batch_comp_star(self, band_pair, auid):
        batches = self.batch_and_sequence_band_pair(band_pair)
        result = batches[batches['auid'] == auid]['batch_id', band_pair[0], band_pair[1],
                                                 self.instr(band_pair[0]), self.instr(band_pair[1])]
        result.sort('batch_id')
        return result

    def batch_target_star(self, band_pair, auid):
        batches = self.batch_band_pair(band_pair)
        result = batches[batches['auid'] == auid]['batch_id', self.instr(band_pair[0]), self.instr(band_pair[1])]
        result.rename_columns([self.instr(band_pair[0]), self.instr(band_pair[1])],
                              [self.targ(band_pair[0]), self.targ(band_pair[1])])
        result.sort('batch_id')
        return result