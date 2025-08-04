import matplotlib as mpl
from cycler import cycler

default_color = '#a08000'

mpl.colormaps.unregister('nightvis')
cdict = {
    'red': [
        (0.0,    0.0,    0.0),
        (0.5,    0.3,    0.3),
        (1.0,    0.6,    0.6),
    ],
    'green': [
        (0.0,    0.0,    0.0),
        (0.5,    0.0,    0.0),
        (1.0,    0.5,    0.5),
    ],
    'blue': [
        (0.0,    0.0,    0.0),
        (1.0,    0.0,    0.0),
    ],
}
mpl.colormaps.register(
    cmap=mpl.colors.LinearSegmentedColormap('nightvis', cdict))

dark_plot_theme = {'lines.linewidth': 1.7,
                   'lines.antialiased': True,
                   'patch.linewidth': 1.0,
                   'patch.facecolor': 'k',
                   'patch.edgecolor': 'k',
                   'patch.antialiased': True,
                   'image.cmap': 'nightvis',
                   'image.origin': 'upper',
                   'font.size': 12.0,
                   'axes.facecolor': '#000000',
                   'axes.edgecolor': default_color,
                   'axes.linewidth': 1.0,
                   'axes.grid': True,
                   'axes.titlesize': 'small',
                   'axes.labelsize': 'small',
                   'axes.labelcolor': default_color,
                   'axes.axisbelow': True,
                   'xtick.major.size': 0,
                   'xtick.minor.size': 0,
                   'xtick.major.pad': 6,
                   'xtick.minor.pad': 6,
                   'xtick.color': default_color,
                   'xtick.direction': 'in',
                   'ytick.major.size': 0,
                   'ytick.minor.size': 0,
                   'ytick.major.pad': 6,
                   'ytick.minor.pad': 6,
                   'ytick.color': default_color,
                   'ytick.direction': 'in',
                   'legend.fancybox': True,
                   'legend.loc': 'best',
                   'figure.figsize': [8, 6],
                   'figure.facecolor': 'k',
                   'figure.edgecolor': 'k',
                   'figure.subplot.hspace': 0.5,
                   'figure.subplot.bottom': 0.15,
                   'savefig.dpi': 72,
                   'savefig.facecolor': 'k',
                   'savefig.edgecolor': 'k',
                   'savefig.transparent': False,
                   'grid.color': default_color,
                   'text.color': default_color,
                   'axes.prop_cycle': cycler('color', ['#a02000', '#800080', '#c00080', '#8000c0', '#808000', '#c08000', '#c04000'])
                   }
