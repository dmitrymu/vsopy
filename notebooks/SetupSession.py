import os
import sys
module_path = os.path.abspath(os.path.join('../src'))
if module_path not in sys.path:
    sys.path.append(module_path)

import astropy.units as u
from astropy.table import unique, vstack
from itertools import dropwhile
import numpy as np
import ipywidgets as widgets
import ccdproc as ccdp
import itertools as it
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch, ImageNormalize, MinMaxInterval
import logging
from matplotlib.patches import Circle
from IPython.display import display

import vso.reduce

def get_image_fov(image):
    preferred = [10, 20, 30, 60, 120, 180] * u.arcmin
    hdr = image.header
    wcs = image.wcs
    p1 = (0, 0)
    p2 = (hdr['naxis1'], hdr['naxis2'])
    c1 = wcs.pixel_to_world(*p1)
    c2 = wcs.pixel_to_world(*p2)
    fov = c1.separation(c2).to(u.arcmin)
    tail = list(dropwhile(lambda x,: x < fov / np.sqrt(2), preferred))
    return tail[0] if len(tail) > 0 else preferred[-1]

def get_image_stars(sd, image):
    object_name = image.header['object']
    objects = ([object_name] if not sd.is_std_field(object_name)
                else [n for n in sd.std_fields['name']
                    if n == object_name or n.startswith(f"{object_name} ")])

    stars = unique(vstack([sd.get_chart(name,
                                       fov=get_image_fov(image))
                          for name in objects]))
    return stars[image.wcs.footprint_contains(stars['radec2000'])]

class ImagePreview:
    def __init__(self, image_dir, stardata, solved_dir) -> None:
        logging.getLogger('astropy').setLevel(logging.ERROR)
        logging.getLogger('root').setLevel(logging.ERROR)
        band_dirs = [d for d in image_dir.iterdir() if d.is_dir()]
        self.files_ = {d.parts[-1]: ccdp.ImageFileCollection(d) for d in band_dirs}
        self.image_ = None
        self.chart_ = None
        self.stardata_ = stardata
        self.band_ = None
        self.solved_dir_ = solved_dir

        bands = list([d.parts[-1] for d in band_dirs])
        self.wgt_play = widgets.Play(value=0, min=0, max=1, interval=1500, description='Play')
        self.wgt_file = widgets.Dropdown(options=[],
                                    description='Image',
                                    )
        self.wgt_file.layout.width='80%'
        self.wgt_band = widgets.Dropdown(options=bands,
                                    description='Band')

        self.wgt_band.observe(self.band_updated, 'value')
        self.wgt_play.observe(self.play_updated, 'value')

    @property
    def image(self):
        return self.image_

    @property
    def chart(self):
        return self.chart_

    def band_updated(self, *args):
        self.band_ = self.wgt_band.value
        self.wgt_file.options = self.files_[self.band_].summary['file']
        self.wgt_file.value = self.wgt_file.options[0]
        self.wgt_play.max = len(self.files_[self.band_].summary)-1
        self.wgt_play.value = 0

    def play_updated(self, *args):
        self.wgt_file.value = self.files_[self.band_].summary['file'][self.wgt_play.value]


    def plot_preview(self, band, file, alpha):
        ifc = self.files_[band]
        path = ifc.location / self.wgt_file.value
        self.image_  = vso.reduce.load_and_solve(path, self.solved_dir_)
        self.chart_ = get_image_stars(self.stardata_, self.image_ )
        wcs = self.image_ .wcs
        interval = MinMaxInterval()
        vmin, vmax = interval.get_limits(self.image_ .data)
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(alpha))
        fig = plt.figure(figsize=(8.0, 8.0))
        ax = plt.subplot(projection=wcs)
        ax.imshow(self.image_.data,  origin='lower', norm=norm , cmap='gray_r')
        for star, n in zip(self.chart_, it.count()):
            r = 10
            c = wcs.world_to_pixel(star['radec2000'])
            ax.add_patch(Circle(c, radius=r, edgecolor='red', facecolor='none', alpha = .5))
            tx, ty = wcs.pixel_shape[0]+60, wcs.pixel_shape[1] - (n+1) * 60
            ax.text(tx, ty, f"{n+1}: {star['auid']}")
            ax.annotate(f"{n+1}", c)
        ax.set_title(f"{self.image_.header['object']}, {self.band_}, "
                     f"UT {self.image_.header['date-obs'].replace('T', ' ')}")
        plt.show()

    def run(self):
        self.band_updated()

        display(self.wgt_play)
        widgets.interact(self.plot_preview,
                        band=self.wgt_band,
                        file=self.wgt_file,
                        alpha=widgets.FloatLogSlider(value=.001, min=-5, max=1, step=.2,
                                                    description="Stretch Factor"))
