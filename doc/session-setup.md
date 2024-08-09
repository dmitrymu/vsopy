# Session Setup

The first step after images are acquired and place into the
[image archive](filesystem.md#image-archive)
is setting up respective session directory inside the
[working directory](filesystem.md#work-directory).  This is
assisted by [SetupSession notebook](../notebooks/SetupSession.ipynb).

The user should change the following values as needed:

```python
OBJ_NAME='SA38'
SESSION_TAG='20240731'
IMAGE_ROOT = '/srv/public/img'
WORK_ROOT = '/srv/public'
```

After that, run all code cells. The notebook provides two preview sections:

1. Image overview, shows a source image and star positions.
1. Photometry preview, show magnified stars fro photometric sequence.
   Stars are overlaid with aperture and annulus. Change radii as needed.
   Disable stars that cannot be reliably measured (e.g. other star inside
   annulus)

Once apertures are chosen, go to the last section and save the settings.
