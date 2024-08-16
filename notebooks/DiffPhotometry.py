import ipywidgets as widgets
import itertools
from astropy.table import QTable
import vso.phot
import matplotlib.pyplot as plt
import numpy as np
import json
import vso.util
from IPython.display import display

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
        return np.sqrt(np.sum(data[band]['err']**2))

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


    def check_err(self, data, band):
        return np.sqrt(np.sum((data[f'check {band}']['mag'] - data[f'check std {band}']['mag'])**2))

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
            dc = data[f'check {band}']['mag'] - data[f'check std {band}']['mag']
            dc_mean = np.mean(dc)
            dc_std = np.std(dc)
            flt = np.abs(dc - dc_mean) < 2*dc_std
            good = data[flt]
            bad = data[~flt]
            ax.errorbar(good['time'], good[f'{band}']['mag'],
                        yerr=good[f'check {band}']['err'], fmt='o', label=band)
            ax.errorbar(bad['time'], bad[f'{band}']['mag'], yerr=bad[f'check {band}']['err'], fmt='o')


        fig = plt.figure(figsize=(8.0, 6.4))
        ax = plt.subplot()
        plot_band(ax, dph, band[0])
        plot_band(ax, dph, band[1])
        ax.invert_yaxis()
        ax.legend()
        plt.show()


    def run(self):
        widgets.interact(self.plot,
                        b=self.wgt_band_,
                        cm=self.wgt_comp_,
                        ck=self.wgt_check_
        )