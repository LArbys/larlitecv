import os,sys,time

s = time.time()
import ROOT
larlitecv_dir = os.environ['LARLITECV_LIBDIR']
for l in [x for x in os.listdir(larlitecv_dir) if x.endswith('.so')]:
    print "load: ",l
    ROOT.gSystem.Load(l)
from ROOT import larlitecv
print "loading libs: ",time.time()-s


s = time.time()
fman_larlite = larlitecv.LarliteFileManager("ex_databnb_larlite.txt")
fman_larlite.initialize()
print "load larlite fman: ",time.time()-s

s = time.time()
fman_larcv = larlitecv.LarcvFileManager("ex_databnb_larcv.txt")
fman_larcv.initialize()
print "load larcv fman: ",time.time()-s


