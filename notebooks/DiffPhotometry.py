import ipywidgets as widgets
import itertools
from astropy.table import QTable
import vso.phot
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import json
import vso.util
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
    prev_settings = vso.util.Settings(settings_path)


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
    def __init__(self, settings, layout) -> None:
        self.phot_ = QTable.read(layout.photometry_file_path)
        self.chart_ = QTable.read(layout.chart_file_path)
        self.provider_ = vso.phot.DataProvider(
            self.phot_, self.chart_
        )
        self.settings_ = settings
        self.bands_ = vso.util.band_pairs(settings.bands)

        self.wgt_band_ = widgets.Dropdown(
            options=[(f"{b[0]} {b[1]}", n) for b, n in zip(self.bands_, itertools.count())],
            description = 'Band'
            )
        self.wgt_band_.observe(self.band_updated, 'value')

        self.wgt_comp_ = widgets.Select(
            options = self.chart_['auid'],
            rows=len(self.chart_),
            description = 'Comp star',
            layout=widgets.Layout(width='auto')
        )
        self.table_comp_ = None
        self.band_updated()
        self.wgt_comp_.observe(self.comp_updated, 'value')

        self.table_check_ = None
        self.wgt_check_ = widgets.Select(
            options = self.chart_['auid'],
            rows=len(self.chart_),
            description = 'Check star',
            layout=widgets.Layout(width='auto')
        )
        self.comp_updated()

    @property
    def settings(self):
        return self.settings_

    def comp_updated(self, *args):
        band = self.bands_[self.wgt_band_.value]
        comp = self.table_comp_[self.wgt_comp_.value]['auid']
        self.table_check_ = self.table_check_err(band, comp)
        self.table_check_.sort(band[0])
        self.wgt_check_.options=[
            (self.format_label(row, band), n)
            for row, n in zip(self.table_check_, itertools.count())]
        self.wgt_check_.rows=len(self.table_check_)

    def format_label(self, row, band):
        return (f"{row['auid']}: σ{band[0]}={row[band[0]]:.3g}, "
                f"σ{band[1]}={row[band[1]]:.3g},  N={row['N']}")

    def band_updated(self, *args):
        band = self.bands_[self.wgt_band_.value]
        self.table_comp_ = self.table_target_err(band)
        self.table_comp_.sort(band[0])
        self.wgt_comp_.options=[
            (self.format_label(row, band), n)
            for row, n in zip(self.table_comp_, itertools.count())]
        self.wgt_comp_.rows=len(self.table_comp_)

    def target_err(self, data, band):
        return np.sqrt(np.sum(data[band]['err']**2) / len(data))

    def band_target_err(self, band, comp):
        xfm = vso.phot.BatchTransformer(band, comp)
        dph = xfm.calculate(self.provider_)
        if len(dph) > 0:
            err_a = self.target_err(dph, band[0])
            err_b = self.target_err(dph, band[1])
            return err_a, err_b, len(dph)
        else:
            return np.nan, np.nan, 0

    def table_target_err(self, band):
        err = [self.band_target_err(band, comp) for comp in self.chart_['auid']]
        zipped = [(comp, e[0], e[1], e[2])
                  for comp, e in zip(self.chart_['auid'], err)
                  if not np.isnan(e[0])]
        return QTable({
            'auid': [z[0] for z in zipped],
            band[0]: [z[1] for z in zipped],
            band[1]: [z[2] for z in zipped],
            'N': [z[3] for z in zipped],
        })


    def check_err_data(self, data, band):
        return data[f'check {band}']['mag'] - data[f'check std {band}']['mag']

    def check_err(self, data, band):
        return np.sqrt(np.sum(self.check_err_data(data, band)**2) / len(data))

    def band_check_err(self, band, comp, check):
        xfm = vso.phot.BatchTransformer(band, comp, check)
        dph = xfm.calculate(self.provider_)
        if len(dph) > 0:
            err_a = self.check_err(dph, band[0])
            err_b = self.check_err(dph, band[1])
            return err_a, err_b, len(dph)
        else:
            return np.nan, np.nan, 0

    def table_check_err(self, band, comp):
        err = [self.band_check_err(band, comp, check)
               for check in self.chart_['auid']
               if check != comp]
        zipped = [(check, e[0], e[1], e[2])
                  for check, e in zip(self.chart_['auid'], err)
                  if not np.isnan(e[0])]
        return QTable({
            'auid': [z[0] for z in zipped],
            band[0]: [z[1] for z in zipped],
            band[1]: [z[2] for z in zipped],
            'N': [z[3] for z in zipped],
        })


    def plot(self, b=0, cm=0, ck=0):
        band = self.bands_[b]
        comp = self.table_comp_[cm]['auid']
        check = self.table_check_[ck]['auid']
        self.settings_.set_comp(band, comp)
        self.settings_.set_check(band, check)
        xfm = vso.phot.BatchTransformer(band, comp, check)
        dph = xfm.calculate(self.provider_)

        def plot_band(ax, data, band):
            ax.errorbar(data['time'], data[f'{band}']['mag'],
                        yerr=data[f'check {band}']['err']/2, fmt='o',
                        label=band, color=BAND_COLOR[band])
            # dc = data[f'check {band}']['mag'] - data[f'check std {band}']['mag']
            # dc_mean = np.mean(dc)
            # dc_std = np.std(dc)
            # flt = np.abs(dc - dc_mean) < 2*dc_std
            # good = data
            # bad = data[~flt]
            # ax.errorbar(good['time'], good[f'{band}']['mag'],
            #             yerr=good[f'check {band}']['err']/2, fmt='o', label=band)
            # ax.errorbar(bad['time'], bad[f'{band}']['mag'], yerr=bad[f'check {band}']['err']/2, fmt='o')

        star = self.chart_.meta.get('star', '???')

        fig = plt.figure(figsize=(10.24, 24.0))
        gs = fig.add_gridspec(4, 1)
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
        ax3.plot(dph['time'], (dph[band[0]]['err']),
                 '.', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        ax3.plot(dph['time'], (dph[band[1]]['err']),
                 '.', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        ax3.set_xlabel('JD')
        # ax3.plot(dph[f'peak {band[0]}'], (dph[band[0]]['err']),
        #          '.', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        # ax3.plot(dph[f'peak {band[1]}'], (dph[band[1]]['err']),
        #          '.', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        # ax3.set_xlabel('Peak')
        ax3.set_ylabel('Target uncertainty')
        ax3.set_title(f'{star} Photometry errors for ${{{band[0]}}}$ and ${{{band[1]}}}$')

        #ax4 = ax3.twinx()
        ax4 = fig.add_subplot(gs[2, 0])
        ax4.plot(dph['time'], (self.check_err_data(dph, band[0])),
                 '+', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        ax4.plot(dph['time'], (self.check_err_data(dph, band[1])),
                 '+', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        ax4.set_ylabel('Check star error')

        ax5 = fig.add_subplot(gs[3, 0])
        ax5.plot(dph['time'], (dph[f'peak {band[0]}']),
                 '.', label=f'Target {band[0]}', color=BAND_COLOR[band[0]])
        ax5.plot(dph['time'], (dph[f'peak {band[1]}']),
                 '.', label=f'Target {band[1]}', color=BAND_COLOR[band[1]])
        ax5.set_xlabel('JD')
        ax5.set_ylabel('Target peak')

        # ax2 = fig.add_subplot(gs[2, 0])
        # ci = dph[f'{band[0]}']['mag'] - dph[f'{band[1]}']['mag']
        # u_ci = np.sqrt(np.square(dph[f'check {band[0]}']['err']) + np.square(dph[f'check {band[1]}']['err']))
        # ax2.errorbar(dph['time'], ci, yerr=u_ci/2, fmt='o', color='grey')
        # ax2.xaxis.set_major_formatter(fmter)
        # ax2.invert_yaxis()
        # ax2.set_xlabel('JD')
        # ax2.set_ylabel('Magnitude')
        # ax2.set_title(f'{star} Color Index ${{{band[0]}-{band[1]}}}$')

        plt.show()


    def run(self):
        widgets.interact(self.plot,
                        b=self.wgt_band_,
                        cm=self.wgt_comp_,
                        ck=self.wgt_check_
        )