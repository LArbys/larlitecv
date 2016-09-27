import os,sys
from time import time

s = time()
import ROOT
print "ROOT: ",time()-s,"secs"

#s = time()
#from ROOT import larlitecv
#print "larlitecv: ",time()-s,"secs"

s = time()
ROOT.gSystem.Load("liblarlitecv")
print "larlitecv base: ",time()-s,"secs"

s = time()
from ROOT import larlitecv
print "from larlitecv: ",time()-s,"secs"

#s = time()
#from larlitecv import larlitecv
#print "data coord: ",time()-s

print "make data co."
dataco = larlitecv.DataCoordinator()

fin_larlite = "ex_databnb_larlite.txt"
fin_larcv   = "ex_databnb_larcv.txt"

flarlite = open( fin_larlite )
lines =  flarlite.readlines()
for l in lines:
    dataco.add_inputfile( l.strip(), "larlite" )

flarcv = open( fin_larcv )
lines =  flarcv.readlines()
for l in lines:
    dataco.add_inputfile( l.strip(), "larcv" )

dataco.initialize()
