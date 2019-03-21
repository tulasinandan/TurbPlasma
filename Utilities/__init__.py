#from os.path import dirname, basename, isfile
#import glob
#modules = glob.glob(dirname(__file__)+"/*.py")
#__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
#[__import__(i) for i in __all__]
from . import (eplots, multigray, multiplt,
        pgeplots, subs, extract_slice)
