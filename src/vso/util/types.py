from collections import namedtuple

ValErr = namedtuple('ValErr', ['val', 'err'])
ValErrDtype = [('val', 'f4'), ('err', 'f4')]

MagErr = namedtuple('ValErr', ['mag', 'err'])
MagErrDtype = [('mag', 'f4'), ('err', 'f4')]
