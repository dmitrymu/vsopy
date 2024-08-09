# Photometry Extraction

Once [session setup](session-setup.md) completed, photometry data have to be
extracted from the source images.  This is done by the command-line script
[src/pipeline/photometry.py](../src/pipeline/photometry.py). Run it as follows:

```{}
photometry.py [-h] -O OBJECT -t TAG -c COMMON_DIR -i IMAGE_DIR

options:
  -h, --help            show this help message and exit
  -O OBJECT, --object OBJECT
                        Object to be worked
  -t TAG, --tag TAG     Tag (date)
  -c COMMON_DIR, --common_dir COMMON_DIR
                        File tree root
  -i IMAGE_DIR, --image_dir IMAGE_DIR
                        Image directory
```

See [file layout](filesystem.md) description for the explanation of parameters.

The script takes the list of stars to be processed from `cetroids.ecsv` and
apertures from `settings.json` in the session directory.

Source images are converted from ADU to electrons and reduced using calibration
frames (auto-selected).  For each reduced image the script measures flux for
every star in the list and calculate its instrumental magnitude. Uncertainty
is propagated from both source and calibration images up to instrumental
magnitudes.

This script may run for many minutes. Image processing is parallelized to
utilize multiprocessor systems.

The result of the script is stored in the session folder as `photometry.ecsv`.
