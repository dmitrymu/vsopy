# Image Acquisition

(TODO: how to generate KStars catalog for a list of AAVSO-supported targets).

Images are acquired by KStars/EKOS. By default it creates the following
directory structure under the given root:

- target name
  - frame type (light, dark, etc)
    - filter name (for light and flat images)
      - image files

After imaging the data should be moved from EKOS directory to image archive.
The archive directory contains subdirectories for every night of observation.

Calibration frames (bias, dark, and flat) are acquired and archived the same way.
