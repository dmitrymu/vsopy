# Simple Differential Photometry

Instrumental magnitudes [extracted](photometry-extraction.md) from source
images can be used for differential photometry.

The script [src/pipeline/diff_simple.py](../src/pipeline/diff_simple.py)
performs simplified differential photometry following [Bruce L. Gary
method](https://www.aavso.org/sites/default/files/Transforms-Sarty.pdf)

```{}
python3 diff_simple.py [-h] -O OBJECT -t TAG -c COMMON_DIR --observer OBSERVER

options:
  -h, --help            show this help message and exit
  -O OBJECT, --object OBJECT
                        Object to be worked
  -t TAG, --tag TAG     Tag (date)
  -c COMMON_DIR, --common_dir COMMON_DIR
                        File tree root
  --observer OBSERVER   AAVSO observer code

```

TODO: parameterize the list of bands

The script uses instrumental magnitudes from `photometry.ecsv` and standard
magnitudes from `chart.ecsv` in the session directory. The result is written
as a file  `report-simple.txt` in AAVSO Extended Format, ready for upload.
