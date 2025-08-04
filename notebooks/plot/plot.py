import astropy.units as u
import matplotlib.pyplot as plt

#from ..phot import phot
from astropy.nddata import Cutout2D
from astropy.visualization import AsinhStretch, ImageNormalize, MinMaxInterval
from contextlib import contextmanager
from matplotlib.patches import Circle

def iterate_grid(items, grid, fig=None):
    for item, n in zip(items, range(len(items))):
        row = n // grid.ncols
        col = n % grid.ncols
        if fig is None:
            yield (item, grid[row, col])
        else:
            yield (item, fig.add_subplot(grid[row, col]))

@contextmanager
def grid_figure(items, ncols=4, width=6.40):
    nrows = (len(items)+1) // ncols + 1
    fig = plt.figure(figsize=(width, (nrows+1) * width / ncols))
    gs = fig.add_gridspec(nrows, ncols)
    try:
        yield fig, iterate_grid(items, gs)
    finally:
#        plt.tight_layout()
        plt.show()


def default_image_normalizer(data):
    interval = MinMaxInterval()
    vmin, vmax = interval.get_limits(data)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch())


def cutout(fig, image, position, size,
                subplot=None,
                wcs = None,
                norm=default_image_normalizer):
    cutout =  Cutout2D(image.data,
                       position,
                       size,
                       wcs=image.wcs if wcs is None else wcs)

    # ax = plt.subplot(projection=cutout.wcs) if subplot is None else fig.add_subplot(subplot, projection=cutout.wcs)
    ax = plt.subplot() if subplot is None else fig.add_subplot(subplot)
    ax.imshow(cutout.data, origin='lower', norm=norm(cutout.data))
    # ax.coords[0].set_axislabel(position.ra.to_string(u.hour))
    # ax.coords[1].set_axislabel("")
    # ax.coords[0].set_ticks([] * u.degree)
    # ax.coords[1].set_ticks([] * u.degree)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    return (ax, cutout)

def cross(ax, center, size, **kwargs):
    wing = size / 2
    ax.plot([center[0] - wing, center[0] + wing], [center[1], center[1]], **kwargs)
    ax.plot([center[0], center[0]], [center[1] - wing, center[1] + wing], **kwargs)

def xcross(ax, center, size, **kwargs):
    wing = size / 2
    ax.plot([center[0] - wing, center[0] + wing], [center[1] - wing, center[1] + wing], **kwargs)
    ax.plot([center[0] - wing, center[0] + wing], [center[1] + wing, center[1] - wing], **kwargs)

def circle(ax, center, radius, **kwargs):
    ax.add_patch(Circle(center, radius=radius, **kwargs))

