import os,sys

import ROOT
larlitecv_dir = os.environ['LARLITECV_LIBDIR']
for l in [x for x in os.listdir(larlitecv_dir) if x.endswith('.so')]:
    print "load: ",l
    ROOT.gSystem.Load(l)
from ROOT import larlitecv


fman = larlitecv.LarliteFileManager("flist.txt")
fman.initialize()
