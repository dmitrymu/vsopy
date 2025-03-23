import astropy.units as u
import ccdproc as ccdp
import subprocess
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from pathlib import Path

def astap_solver(file_path, solved_dir, radius=10*u.deg):
    file_name = Path(file_path).name
    solver_path = Path(solved_dir) / file_name
    wcs_path = solver_path.with_suffix('.wcs')
    if not wcs_path.exists():
        solver = f"astap_cli -f {file_path} -wcs -sip -r {radius.value} -o {solver_path}"
        rc = subprocess.run(solver, shell=True )
        if rc.returncode != 0:
            raise RuntimeError(f"ASTAP solver failed for {file_path}")
    hdul_wcs = fits.open(wcs_path)
    return hdul_wcs[0].header

def update_wcs(image, wcs_header):
    header = fits.Header(image.header)
    header.update(wcs_header)
    image.wcs = WCS(header)
    return image

def load_and_solve(file_path, solved_dir, radius=10*u.deg):
    file_name = Path(file_path).name
    solver_path = Path(solved_dir) / file_name
    wcs_path = solver_path.with_suffix('.wcs')
    if not wcs_path.exists():
        solver = f"astap_cli -f {file_path} -wcs -sip -r {radius.value} -o {solver_path}"
        rc = subprocess.run(solver, shell=True )
        if rc.returncode != 0:
            raise RuntimeError(f"ASTAP solver failed for {file_path}")
    hdul_wcs = fits.open(wcs_path)
    image = CCDData.read(file_path, unit='adu')
    header = fits.Header(image.header)
    header.update(hdul_wcs[0].header)
    image.wcs = WCS(header)
    return image
