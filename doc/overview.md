# VSO.py Framework

The goal of the framework is bulk processing of photometry data and
visual tools to control the data quality.

## Dependencies

The framework is designed to work with images acquired by the open-source
software [KStars](https://kstars.kde.org/) . VSO.py relies on directory
structure and FITS keywords produced by KStars (namely its EKOS subsystem).

VSO.py uses reference photomery data from [AAVSO](https://www.aavso.org/)
and produces observation reports in
[AAVSO Extended File Format](https://www.aavso.org/aavso-extended-file-format),
ready for uploading to [WevbObs](https://www.aavso.org/webobs).

The framework uses command-line version of [ASTAP](https://www.hnsky.org/astap.htm)
as a local plate solver.

## Astropy Features

VSO.py is based on [astropy](https://docs.astropy.org/en/stable/index.html) and
affiliated libraries like [ccdproc](https://ccdproc.readthedocs.io/en/latest/index.html)
and [photutils](https://photutils.readthedocs.io/en/stable/index.html).

One astropy feature of particular importance is thorough error propagation during
image processing.

## Visual Interface

VSO.py includes some [Jupyter notebooks](https://jupyter.org/) as visual
interface for configuration and data assessment. [Matplotlib](https://matplotlib.org/)
is used for data visualization and [ipywidgets](https://ipywidgets.readthedocs.io/)
provide some limited level of interactivity.

[Visual Studio Code](https://code.visualstudio.com/) is the recommended interface
for running VSO.py.
[Jupyter extension](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter)
is required to work with notebooks.
