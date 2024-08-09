# VSO.py

This is a loose set of tools for processing variable star photometry data.
It is based on [Astropy](https://www.astropy.org/) and affiliated packages.
Jupyter notebooks are used as a visual interface. Batch processing is
implemented as a set command-line utilities.  This is rather research project
then released software.

See [documentation](doc/_toc.md) for details.

## Functionality

1. Providing file system layout for image archive and working data.
1. Creating master calibration images from original bias, dark, and flat
   frames.
1. Applying calibration to light frames.
1. Plate solving (currently using [ASTAP](https://www.hnsky.org/astap.htm)
   solver).
1. Interface to [AAVSO](https://www.aavso.org) photometry databases.
1. Aperture photometry.
1. Differential photometry (simplified)
1. Visual assessment of images.

## Installation

Pre-requisites: git, Python 3, venv.

Clone the repository from [GitHub](https://github.com/dmitrymu/vsopy).

Initialize repository: from the repository directory run `./init.sh`.

Alternatively from VS code:
`Ctrl+Shift+P` -> `Python: Create Environment...` and mark checkbox for
`requirements.txt` when prompted to select dependencies. Then run
`git config --local core.hooksPath .githooks/` from terminal
