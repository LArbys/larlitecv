import ROOT,os
if not 'LARLITECV_BASEDIR' in os.environ:
    print '$LARLITECV_BASEDIR shell env. var. not found (run configure.sh)'
    raise ImportError
#force loading larlite libs FIRST because auto-loading in ROOT6 does not properly work
if 'LARLITE_BASEDIR' in os.environ:
    from larlite import larlite
if 'LARCV_BASEDIR' in os.environ:
    from larcv import larcv
larlitecv_dir = os.environ['LARLITECV_LIBDIR']
for l in [x for x in os.listdir(larlitecv_dir) if x.endswith('.so')]:
    ROOT.gSystem.Load(l)
from ROOT import flashana
flashana.LightPath()
from ROOT import larlitecv


