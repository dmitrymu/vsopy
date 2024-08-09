import astropy.units as u
import ccdproc as ccdp
import numpy as np
import logging
import tempfile
import psutil
from astropy.nddata import CCDData
from astropy.stats import mad_std
from pathlib import Path
from vso.reduce import CalibrationMatcher
from vso.data import CameraRegistry
from vso.util import FrameType

def round_Mb(x):
    return (x >> 20) << 20

class MasterBuilder:
    def __init__(self, type, image_dir, output_dir, tmp_dir, overwrite=False, delete_tmp=True) -> None:
        self.type_ = type
        self.image_dir_ = Path(image_dir)
        self.output_dir_ = Path(output_dir)
        self.tmp_dir_ = Path(tmp_dir)
        self.overwrite_ = overwrite
        self.delete_tmp_ = delete_tmp

        self.ifc = [ccdp.ImageFileCollection(filter_dir)
                    for filter_dir in self.image_dir_.iterdir() if filter_dir.is_dir()] \
            if self.type_ == FrameType.FLAT \
            else [ccdp.ImageFileCollection(self.image_dir_)]

        self.matcher = CalibrationMatcher(self.output_dir_)
        logging.getLogger('astropy').setLevel(logging.ERROR)
        logging.getLogger('root').setLevel(logging.ERROR)

    def deviate_correct(image, gain=None, readnoise=None):
        return ccdp.gain_correct(ccdp.create_deviation(image,
                                    gain,
                                    readnoise), gain)

    def process_image(self, path, tmp_dir):
        image = CCDData.read(path, unit='adu')
        camera_name = image.header['instrume']
        camera = CameraRegistry.get(camera_name)
        image.divide(camera.adu_scale)
        image_gain = image.header['gain']
        e_gain = camera.gain_to_e(image_gain)
        e_noise = camera.read_noise(image_gain)
        reduced = ccdp.gain_correct(
            ccdp.create_deviation(image,
                                  e_gain,
                                  e_noise),
            e_gain)
        cal = self.matcher.match(image.header)

        if not cal.bias:
            reduced.meta['bias-sub'] = 'F'
        else:
            reduced = ccdp.subtract_bias(reduced, cal.bias)
            reduced.meta['bias-sub'] = 'T'

        if not cal.dark:
            reduced.meta['dark-sub'] = 'F'
        else:
            reduced = ccdp.subtract_dark(reduced, cal.dark,
                                         exposure_time='exptime',
                                         exposure_unit=u.second,
                                         scale=True)
            reduced.meta['dark-sub'] = 'T'
        reduced.write(Path(tmp_dir) / path.parts[-1])


    def prepare_images(self, group, tmp_dir):
        paths= [self.image_dir_ / r['file'] for r in group]
        for p in paths:
            self.process_image(p, tmp_dir)

    def create_master(self, num_frames, keys, temp, tmp_dir):
        deviated = ccdp.ImageFileCollection(tmp_dir)
        mem_limit = round_Mb((psutil.virtual_memory().available*2)//3)
        print(f"using {mem_limit/1024/1024} MB of RAM")
        master = ccdp.combine(deviated.files_filtered(include_path=True),
                              method='average',
                              scale=lambda x: 1 /
                              np.median(
                                  x) if self.type_ == FrameType.FLAT else 1,
                              sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                              mem_limit=mem_limit)
        master.meta['combined'] = 'T'
        master.meta['frame-ct'] = num_frames
        # if 'darktime' in master.meta:
        #     master.meta['darktime'] = master.meta['darktime'] * master.meta['frame-ct']
        master.meta['comment'] = 'Created by VSO master image pipeline'
        exp_descr = '' if self.type_ != FrameType.DARK else f"_e{keys['exptime']}"
        result_name = (f"master_{self.type_}"
                        # f"{'-b' if self.bias else ''}"
                        # f"{'-d' if self.dark else ''}"
                        f"{'_' + keys['filter'] if self.type_ == FrameType.FLAT else ''}"
                        f"_g{keys['gain']:g}"
                        f"{exp_descr}"
                        f"_o{keys['offset']}"
                        f"_b{keys['xbinning']}x{keys['ybinning']}"
                        f"_t{temp:.3g}"
                        # f"_{self.args.tag}"
                        ".fits")
        master.write(self.output_dir_ / result_name,
                     overwrite=self.overwrite_)

    def process(self):

        columns = ['instrume', 'gain', 'xbinning', 'ybinning', 'offset']
        if self.type_ == FrameType.DARK:
            columns.append('exptime')
        if self.type_ == FrameType.FLAT:
            columns.append('filter')

        for ifc in self.ifc:
            grouped = ifc.summary.group_by(columns).groups
            for group, keys in zip(grouped, grouped.keys):
                print(dict(keys))
                with tempfile.TemporaryDirectory(dir=self.tmp_dir_,
                                                 delete=self.delete_tmp_) as tmp:
                    try:
                        tmp_dir = self.tmp_dir_/tmp
                        self.prepare_images(group, tmp_dir)
                        self.create_master(len(group), keys, np.mean(
                            group['ccd-temp']), tmp_dir)
                    except Exception as e:
                        print(f"\nFailed: {e}")
