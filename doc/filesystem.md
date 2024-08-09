# File System Layout

VSO.py enforces some conventions on how input and output data are placed
in the file system. There are two separate layouts, the image archive directory
and the work directory.

## Image Archive

This a directory tree holding original images. The structure is:

```{}
<root>/<tag>/<target>/<frame type>/<filter>
```

for example: `/srv/public/20230122/RR_Lyr/Light/V`

Here `<tag>` may be, in fact, some arbitrary part of path, e.g. `2023/01/22`
or `current`.  Typically it holds observation date. Directory structure
below `<tag>` level was chosen to fit [KStars](acquisition.md) default structure.

## Work Directory

This is a separate directory containing data produced and consumed by various
parts of VSO.py. Its top-level subdirectories are:

- `calibr`, stores [master calibration frames](calibration.md)
- `charts`, stores [reference data](charts.md). e.g. photometry from AAVSO.
- `session`, stores processing results, see [below](#session-directory)
- `tmp`, stores temporary files

## Session Directory

Session is a logical unit of work, e.g. all observations of a given target
during single night. It corresponds to the part of image archive under
`<root> / <tag> / <target>` and has the same structure:

```{}
<work directory>/session/<tag>/<target>
```

For example, results of processing images from `/srv/public/20230122/RR_Lyr`
may be placed under `/home/user/work/session/20230122/RR_Lyr`.
