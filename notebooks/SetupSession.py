import os
import sys
module_path = os.path.abspath(os.path.join('../src'))
if module_path not in sys.path:
    sys.path.append(module_path)

import astropy.units as u
from astropy.table import unique, vstack, QTable
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
import json
import vso.reduce
import vso.phot
import vso.plot as vp
import matplotlib.pyplot as plt

def use_dark_theme(**kwargs):
    theme = dict(vp.dark_plot_theme)
    theme.update(kwargs)
    plt.style.use(theme)

def confirm_settings(preview, layout):
        # settings, settings_path):
    settings = preview.settings
    settings_path = layout.settings_file_path
    prev_settings = {}
    if settings_path.exists():
        with open(settings_path) as file:
            prev_settings = json.load(file)


    button = widgets.Button(description='Save settings')
    text_layout = widgets.Layout(width='50%',
                                height='16em')
    prev = widgets.Textarea(value=json.dumps(prev_settings, indent=2),
                            description='Previous',
                            layout=text_layout,
                            disabled=True)
    curr = widgets.Textarea(value=json.dumps(settings, indent=2),
                            description='Current',
                            layout=text_layout,
                            disabled=True)
    done = widgets.Valid(value=False, description='Not saved', readout='')
    def on_button_clicked(b):
        with open(settings_path, mode='w') as file:
            json.dump(settings, file)
        preview.centroid.write(layout.centroid_file_path, format='ascii.ecsv', overwrite=True)
        done.value=True
        done.description='Saved'

    button.on_click(on_button_clicked)
    box_layout=widgets.Layout(width='90%',
                                height='auto')
    display(widgets.VBox([widgets.HBox([prev, curr]), widgets.HBox([button, done])],
                         layout=box_layout))
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
    def __init__(self, image_dir, stardata, blacklist_path, solved_dir) -> None:
        logging.getLogger('astropy').setLevel(logging.ERROR)
        logging.getLogger('root').setLevel(logging.ERROR)
        band_dirs = [d for d in image_dir.iterdir() if d.is_dir()]
        self.files_ = {d.parts[-1]: ccdp.ImageFileCollection(d) for d in band_dirs}
        self.image_ = None
        self.chart_ = None
        self.stardata_ = stardata
        self.band_ = None
        self.solved_dir_ = solved_dir
        self.blacklist_path = blacklist_path
        if self.blacklist_path.exists():
            with open(self.blacklist_path) as file:
                blacklist = json.load(file)
                self.blacklist = set(blacklist['blacklist']) if 'blacklist' in blacklist else []
        else:
            self.blacklist=set()

        bands = list([d.parts[-1] for d in band_dirs])
        self.wgt_play = widgets.Play(value=0, min=0, max=1, interval=500, description='Play')
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

    def highlight_star(self, ax, pos, title, r ):
        ax.add_patch(Circle(pos, radius=r, edgecolor='orange', facecolor='none', alpha = .5))
        # tx, ty = wcs.pixel_shape[0]+60, wcs.pixel_shape[1] - (n+1) * 60
        # ax.text(tx, ty, f"{n+1}: {star['auid']}")
        # ax.annotate(f"{n+1}", c)

    def plot_preview(self, band, file, alpha):
        ifc = self.files_[band]
        path = ifc.location / file
        if str(path) in self.blacklist:
            display(widgets.Label(f"Blacklisted file {path}"))
            return
        try:
            self.image_  = vso.reduce.load_and_solve(path, self.solved_dir_)
            self.chart_ = get_image_stars(self.stardata_, self.image_ )
            wcs = self.image_ .wcs
            interval = MinMaxInterval()
            vmin, vmax = interval.get_limits(self.image_ .data)
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(alpha))
            fig = plt.figure(figsize=(10.24, 10.24))
            #ax = plt.subplot(projection=wcs)
            ax = plt.subplot()
            ax.imshow(self.image_.data,  origin='lower', norm=norm)
            for star, n in zip(self.chart_, it.count()):
                r = 20
                c = wcs.world_to_pixel(star['radec2000'])
                self.highlight_star(ax, c, f"{n+1}", r)
            ax.set_title(f"{self.image_.header['object']}, {self.band_}, "
                        f"UT {self.image_.header['date-obs'].replace('T', ' ')}")
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.show()
        except:
            self.blacklist.add(str(path))
            with open(self.blacklist_path, mode='w') as file:
                json.dump(dict(blacklist=list(self.blacklist)), file)

    def run(self):
        self.band_updated()

        display(self.wgt_play)
        widgets.interact(self.plot_preview,
                        band=self.wgt_band,
                        file=self.wgt_file,
                        alpha=widgets.FloatLogSlider(value=.001, min=-5, max=1, step=.2,
                                                    description="Stretch Factor"))

class PreviewPhotometry:
    def __init__(self, image, chart, settings_path, centroid_path) -> None:
        self.r_ap = 5*u.arcsec
        self.r_ann = (10*u.arcsec, 15*u.arcsec)
        self.stars_disabled = set()
        self.settings_ = {}
        if settings_path.exists():
            with open(settings_path) as file:
                self.settings_ = json.load(file)
                if 'aperture' in self.settings_:
                    a = self.settings_['aperture']
                    unit = u.Unit(a['unit'])
                    self.r_ap = a['r_ap'] * unit
                    self.r_ann = a['r_in'] * unit, a['r_out'] * unit
                if 'disabled' in self.settings_:
                    self.stars_disabled = set(self.settings_['disabled'])

        self.image_ = image.divide(4)
        self.image_.header = image.header
        draft = chart['auid', 'radec2000']
        draft.rename_column('radec2000', 'sky_centroid')
        if centroid_path.exists():
            self.centroid = QTable.read(centroid_path, format='ascii.ecsv')
        else:
            self.centroid = vso.phot.measure_photometry(image, draft, self.r_ap, self.r_ann)['auid', 'sky_centroid']

        self.wgt_enable = [widgets.Checkbox(value=star['auid'] not in self.stars_disabled,
                                            description=star['auid'],
                                            indent=False,
                                            layout=widgets.Layout(width='auto')
                                            ) for star in chart]
        def cb_updated(event):
            owner = event['owner']
            star = owner.description
            if owner.value:
                self.stars_disabled.remove(star)
            else:
                self.stars_disabled.add(star)

        for wgt in self.wgt_enable:
            wgt.observe(cb_updated, 'value')


    @property
    def pixel_scale(self):
        return self.image_.header['SCALE'] * u.arcsec / u.pix

    @property
    def settings(self):
        result=self.settings_
        result['aperture']['unit'] = str(u.arcsec)
        result['aperture']['r_ap'] = float(self.r_ap.value)
        result['aperture']['r_in'] = float(self.r_ann[0].value)
        result['aperture']['r_out'] = float(self.r_ann[1].value)
        result['disabled'] = [str(s) for s in self.stars_disabled]
        # return {
        #     "aperture": {
        #         "unit": str(u.arcsec),
        #         "r_ap": float(self.r_ap.value),
        #         "r_in": float(self.r_ann[0].value),
        #         "r_out": float(self.r_ann[1].value)
        #     },
        #     "disabled": [str(s) for s in self.stars_disabled]
        # }
        return result

    def stars_photometry(self, r, r_in=None, r_out=None, tile=60, ncols=6, width=12.80):
        self.r_ap = r*u.arcsec
        self.r_ann = (r_in*u.arcsec, r_out*u.arcsec)
        self.ph = vso.phot.measure_photometry(self.image_, self.centroid, self.r_ap, self.r_ann)
        nrows = (len(self.ph)+1) // ncols + 1
        for row in range(nrows):
            sl = slice(row*ncols, (row+1)*ncols)
            wgt_grid = widgets.GridBox(self.wgt_enable[sl], layout=widgets.Layout(
                grid_template_columns=f"repeat({ncols}, {1/ncols:.2%})", width=f'min(100%, {int(width*100)}px)'))
            display(wgt_grid)
            with vp.grid_figure(self.ph[sl], ncols=ncols, width=width) as (fig, cells):
                #fig.suptitle(f"Stars for object {self.image_.header['OBJECT']}")
                for star, cell in cells:
                    ax, cut = vp.cutout(fig, self.image_, star['sky_centroid'], (tile, tile),
                                    subplot=cell)
                    ax.set_title(f"{star['name'] if 'name' in star.keys() else star['auid']}\n"
#                                fr"$\ {self.image_.header['FILTER']}_{{instr}} = {star['M']['mag'].value:.3g}$"
                                '\n'
                                f"SNR {star['snr']:.3g}\n"
                                f"Peak {star['peak']:.2%}")
                    if star['auid'] in self.stars_disabled:
                        vp.xcross(ax, cut.center_cutout, tile-2,
                                color='orange', alpha=.5)
                    c = cut.wcs.world_to_pixel(star['sky_centroid'])
                    r_ap = (r/self.pixel_scale).value
                    r_1 = (r_in/self.pixel_scale).value if r_in is not None else r_ap * 2
                    r_2 = (r_out/self.pixel_scale).value if r_out is not None else r_ap * 3
                    vp.circle(ax, c, r_ap,
                              linewidth=2, edgecolor='orange', facecolor='none', alpha=.2)
                    vp.circle(ax, c, r_1,
                              linewidth=2, edgecolor='orange', facecolor='none', alpha=.2)
                    vp.circle(ax, c, r_2,
                              linewidth=2, edgecolor='orange', facecolor='none', alpha=.2)

    def run(self, width=12.80):
        widgets.interact(self.stars_photometry,
                        r=widgets.FloatText(value=self.r_ap.value, description='aperture "'),
                        r_in=widgets.FloatText(value=self.r_ann[0].value, description='annulus inner "'),
                        r_out=widgets.FloatText(value=self.r_ann[1].value, description='annulus outer "'),
                        tile=widgets.IntSlider(value=40, min=20, max=60, description="Tile size"),
                        ncols=widgets.IntSlider(value=6, min=1, max=12, description="TNumber of columns"),
                        width=width
                        )
