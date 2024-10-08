# Calibration Frames

Bias, dark, and light frames should be combined into master frames.

There is [`src/pipeline/create_master.py`](../src/pipeline/create_master.py)
command-line script for this purpose.  Use as follows:

```{}
python3 create_master.py [-h] -w WORK_DIR [-i IMAGE_DIR] [-m MEMORY_LIMIT] [--no-cleanup] [--overwrite]

Traverse image directory and create master frames from all discovered calibration images

options:
  -h, --help            show this help message and exit
  -w WORK_DIR, --work-dir WORK_DIR
                        Work directory
  -i IMAGE_DIR, --image-dir IMAGE_DIR
                        Image directory
  -m MEMORY_LIMIT, --memory-limit MEMORY_LIMIT
                        memory limit in MB
  --no-cleanup          Do not remove temporary files
  --overwrite           Overwrite output files
```

Here `IMAGE_DIR` is a directory containing calibration frames, e.g.
`<root>/20240120/Calibr/Dark`.  `COMMON_DIR` is a root working directory.
Master frames are written to `COMMON_DIR/calibr`.  `COMMON_DIR/tmp` is
used to hold temporary files.

Calibration frames are much larger than source images.  The reasons are:

- they are stored in float format because of ADU to electrons conversion;
- they contain per-pixel uncertainty data.

Combiner algorithm is memory-hungry, here comes the need for temporary files
(see [ccdproc.combine](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.combine.html#ccdproc.combine)
for details).
