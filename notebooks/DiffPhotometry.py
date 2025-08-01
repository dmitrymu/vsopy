import ipywidgets as widgets
import itertools
from astropy.table import QTable, join
import vsopy.phot
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import json
import vsopy.util
import astropy.units as u
from IPython.display import display

BAND_COLOR = dict(B='blue',
                  V='green',
                  Rc='red',
                  Ic='magenta')

class ScalarFormatter(ticker.ScalarFormatter):
    def __init__(self, useOffset=None, useMathText=None, useLocale=None):
        ticker.ScalarFormatter.__init__(self, useOffset, useMathText, useLocale)


    def get_offset(self, txt=''):
        if self.orderOfMagnitude:
            txt += u'\u00D7' + '1e%d' % self.orderOfMagnitude + ' '

        if self.offset:
            if self.offset > 0:
                txt += '+'
            txt += f'{self.offset:,.0f}'  # self.format_data(self.offset)

        return self.fix_minus(txt)

def confirm_settings(preview, layout):
        # settings, settings_path):
    settings = preview.settings
    settings_path = layout.settings_file_path
    prev_settings = vsopy.util.Settings(settings_path)


    button = widgets.Button(description='Save settings')
    text_layout = widgets.Layout(width='50%',
                                height='32em')
    prev = widgets.Textarea(value=json.dumps(prev_settings.data, indent=2),
                            description='Previous',
                            layout=text_layout,
                            disabled=True)
    curr = widgets.Textarea(value=json.dumps(settings.data, indent=2),
                            description='Current',
                            layout=text_layout,
                            disabled=True)
    done = widgets.Valid(value=False, description='Not saved', readout='')
    def on_button_clicked(b):
        with open(settings_path, mode='w') as file:
            json.dump(settings.data, file)
        done.value=True
        done.description='Saved'

    button.on_click(on_button_clicked)
    box_layout=widgets.Layout(width='90%',
                                height='auto')
    display(widgets.VBox([widgets.HBox([prev, curr]), widgets.HBox([button, done])],
                         layout=box_layout))

class PreviewDiffPhotometry:
    def __init__(self, settings, layout, save_tables=False) -> None:
        self.save_tables_ = save_tables
        self.session_path = layout.root_dir
        self.provider_ = vsopy.phot.BatchDataProvider(layout)
        self.settings_ = settings
        self.bands_ = vsopy.util.band_pairs(settings.bands)
        self.xfm_ = None

        self.wgt_band_ = widgets.Dropdown(
            options=[(f"{b[0]} {b[1]}", n) for b, n in zip(self.bands_, itertools.count())],
            description = 'Band'
            )
        self.wgt_band_.observe(self.band_updated, 'value')

        self.auids_ = []
        self.wgt_comp_ = widgets.Select(
            options = self.auids_,
            rows=len(self.auids_),
            layout=widgets.Layout(width='95%')
        )
        self.table_comp_ = None
        self.wgt_comp_.observe(self.comp_updated, 'value')

        self.table_check_ = None
        self.wgt_check_ = widgets.Select(
            options = self.auids_,
            rows=len(self.auids_),
            layout=widgets.Layout(width='95%')
        )
        self.band_updated()

        imin = 0
        imax = len(self.provider_.batches_)-1
        self.wgt_time_ = widgets.IntRangeSlider(min=imin, max=imax, value=[imin, imax],
                                                description = 'Time range',
                                                readout=False,
                                                layout=widgets.Layout(width='50%'))
        self.wgt_range_ = widgets.Label("time", layout=widgets.Layout(width='50%'))
        self.wgt_time_.observe(self.time_updated, 'value')
        self.time_updated()
        self.wgt_ui_ = widgets.VBox([self.wgt_band_,
                                      widgets.HBox([widgets.VBox([widgets.Label("Comp star"), self.wgt_comp_], layout=widgets.Layout(width='50%')),
                                                    widgets.VBox([widgets.Label("Check star"), self.wgt_check_], layout=widgets.Layout(width='50%'))],
                                                    layout=widgets.Layout(width='100%')),
                                      widgets.HBox([self.wgt_time_, self.wgt_range_])])

    @property
    def settings(self):
        return self.settings_

    def time_updated(self, *args):
        imin, imax = self.wgt_time_.value
        self.wgt_range_.value = f"UTC {self.provider_.batches_['time'][imin]} - {self.provider_.batches_['time'][imax]}"
        band = self.bands_[self.wgt_band_.value]
        self.settings_.photometry(band).set_start(self.provider_.batches_['time'][imin].jd)
        self.settings_.photometry(band).set_finish(self.provider_.batches_['time'][imax].jd)

    def comp_updated(self, *args):
        band = self.bands_[self.wgt_band_.value]
        comp = self.table_comp_[self.wgt_comp_.value]['auid']
        self.table_check_ = self.table_check_err(band, comp)
        self.table_check_.sort(band[0])
        self.wgt_check_.options=[
            (self.format_label(row, band), n)
            for row, n in zip(self.table_check_, itertools.count())]
        self.wgt_check_.value = 0
        self.wgt_check_.rows=len(self.table_check_)

    def format_label(self, row, band):
        return (f"{row['auid']}: Σ{band[0]}={row[band[0]]:.3g}, "
                f"Σ{band[1]}={row[band[1]]:.3g},  N={row['N']}")

    def band_updated(self, *args):
        band = self.bands_[self.wgt_band_.value]
        self.xfm_ = vsopy.phot.batch_create_simple_transform(self.provider_, band)
        self.auids_ = self.provider_.sequence_band_pair(band)['auid']
        self.table_comp_ = self.table_target_err(band)
        self.table_comp_.sort(band[0])
        self.wgt_comp_.options=[
            (self.format_label(row, band), n)
            for row, n in zip(self.table_comp_, itertools.count())]
        self.wgt_comp_.rows=len(self.table_comp_)
        self.wgt_comp_.value = 0
        self.comp_updated()

    def target_err(self, data, band):
        return np.sqrt(np.sum(data[band]['err']**2) / len(data))

    def band_target_err(self, band, comp):
        dph = vsopy.phot.batch_apply_simple_transform(self.provider_, self.xfm_, band, comp, self.provider_.target_auid)
        if len(dph) > 0:
            err_a = self.target_err(dph, band[0])
            err_b = self.target_err(dph, band[1])
            return err_a, err_b, len(dph)
        else:
            return np.nan, np.nan, 0

    def table_target_err(self, band):
        err = [self.band_target_err(band, comp) for comp in self.auids_]
        zipped = [(comp, e[0], e[1], e[2])
                  for comp, e in zip(self.auids_, err)
                  if not np.isnan(e[0])]
        return QTable({
            'auid': [z[0] for z in zipped],
            band[0]: [z[1] for z in zipped],
            band[1]: [z[2] for z in zipped],
            'N': [z[3] for z in zipped],
        })


    def check_err_data(self, data, band):
        return data[f'check {band}']['mag'] - data[f'check std {band}']['mag']

    def check_err(self, data, goal):
        return np.sqrt(np.sum((data - goal)**2) / len(data))

    def band_check_err(self, band, comp, check):
        dph = vsopy.phot.batch_apply_simple_transform(self.provider_, self.xfm_, band, comp, check)
        check_data = self.provider_.check_band_pair(band, check)
        if len(dph) > 0:
            err_a = self.check_err(dph[band[0]]['mag'], check_data[band[0]]['mag'])
            err_b = self.check_err(dph[band[1]]['mag'], check_data[band[1]]['mag'])
            return err_a, err_b, len(dph)
        else:
            return np.nan, np.nan, 0

    def table_check_err(self, band, comp):
        auids = self.auids_.value.tolist()
        auids.remove(comp)
        err = [self.band_check_err(band, comp, check)
               for check in auids]
        zipped = [(check, e[0], e[1], e[2])
                  for check, e in zip(auids, err)
                  if not np.isnan(e[0])]
        return QTable({
            'auid': [z[0] for z in zipped],
            band[0]: [z[1] for z in zipped],
            band[1]: [z[2] for z in zipped],
            'N': [z[3] for z in zipped],
        })


    def plot(self, b=0, cm=0, ck=0, tr=(0,1)):
        band = self.bands_[b]
        comp = self.table_comp_[cm]['auid']
        check = self.table_check_[ck]['auid']
        check_data = self.provider_.check_band_pair(band, check)

        self.settings_.photometry(band).set_comp(comp)
        self.settings_.photometry(band).set_check(check)
        self.settings_.photometry(band).set_start(self.provider_.batches_['time'][tr[0]].jd)
        self.settings_.photometry(band).set_finish(self.provider_.batches_['time'][tr[1]].jd)
        batches = self.provider_.batches_[tr[0]:tr[1]+1]
        dph = join(vsopy.phot.batch_apply_simple_transform(self.provider_, self.xfm_, band, comp), batches, 'batch_id')
        dph_check = join(vsopy.phot.batch_apply_simple_transform(self.provider_, self.xfm_, band, comp, check), batches, 'batch_id')

        def plot_band(ax, data, band):
            ax.errorbar(data['time'].jd, data[f'{band}']['mag'],
                        yerr=data[f'{band}']['err']/2, fmt='.',
                        label=band, color=BAND_COLOR[band])

        star = self.provider_.sequence_.meta.get('star', '???')

        fig = plt.figure(figsize=(10.24, 10.24))
        gs = fig.add_gridspec(3, 1)
        #ax = plt.subplot()
        #flt = dph[f'check {band[0]}']['err'] < 1 * u.mag
        ax = fig.add_subplot(gs[0, 0])
        plot_band(ax, dph, band[0])
        plot_band(ax, dph, band[1])
        ax.invert_yaxis()
        fmter = ScalarFormatter()
        fmter.set_powerlimits((-4, 8))
        ax.xaxis.set_major_formatter(fmter)
        ax.legend()
        ax.set_xlabel('JD')
        ax.set_ylabel('Magnitude')
        ax.set_title(f'{star} Light Curves ${{{band[0]}}}$ and ${{{band[1]}}}$')

        ax3 = fig.add_subplot(gs[1, 0])
        ax3.xaxis.set_major_formatter(fmter)
        ax3.plot(dph['time'].jd, (dph[band[0]]['err']),
                 '.', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        ax3.plot(dph['time'].jd, (dph[band[1]]['err']),
                 '.', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        ax3.set_xlabel('JD')
        ax3.set_ylabel('Target uncertainty')
        ax3.set_title(f'{star} Photometry errors for ${{{band[0]}}}$ and ${{{band[1]}}}$')

        ax4 = fig.add_subplot(gs[2, 0])
        ax4.plot(dph['time'].jd, ((dph_check[band[0]]['mag'] - check_data[band[0]]['mag'])),
                 '+', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        ax4.plot(dph['time'].jd, ((dph_check[band[1]]['mag'] - check_data[band[1]]['mag'])),
                 '+', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        ax4.set_ylabel('Check star error')

        plt.tight_layout()
        plt.show()


    def run(self):
        w = widgets.interactive_output(self.plot, dict(
                        b=self.wgt_band_,
                        cm=self.wgt_comp_,
                        ck=self.wgt_check_,
                        tr=self.wgt_time_)
        )
        display(self.wgt_ui_, w)