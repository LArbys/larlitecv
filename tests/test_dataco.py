import os,sys
from time import time

s = time()
import ROOT
from ROOT import std
print "ROOT: ",time()-s,"secs"

#s = time()
#from ROOT import larlitecv
#print "larlitecv: ",time()-s,"secs"

#s = time()
#ROOT.gSystem.Load("liblarlitecv")
#print "larlitecv base: ",time()-s,"secs"

#s = time()
#from ROOT import larlitecv
#print "from larlitecv: ",time()-s,"secs"

#s = time()
#from larlitecv import larlitecv
#print "data coord: ",time()-s
s = time()
from larlite import larlite
from larcv import larcv
from larlitecv import larlitecv
print "load larlitecv: ",time()-s

print "make data co."
dataco = larlitecv.DataCoordinator()

fin_larlite = "ex_databnb_larlite.txt"
fin_larcv   = "ex_databnb_larcv.txt"

flarlite = open( fin_larlite )
lines =  flarlite.readlines()
for l in lines:
    dataco.add_inputfile( l, "larlite" )

flarcv = open( fin_larcv )
lines =  flarcv.readlines()
for l in lines:
    dataco.add_inputfile( l, "larcv" )

dataco.configure( "datacoord.cfg", "StorageManager", "IOManager" )

dataco.initialize()

print "entries (larlite): ",dataco.get_nentries("larlite")
print "entries (larcv):   ",dataco.get_nentries("larcv")

larlite_io = dataco.get_larlite_io()
larcv_io   = dataco.get_larcv_io()

for i in range(0,dataco.get_nentries("larlite")):
    dataco.goto_entry( i, "larlite" )
    event_imgs = larcv_io.get_data( larcv.kProductImage2D, "tpc" ) 
    print "Entry: ",i
    print " LARLITE:",larlite_io.run_id(),larlite_io.subrun_id(),larlite_io.event_id()
    print " LARCV: ",event_imgs.run(),event_imgs.subrun(),event_imgs.event() 

    
    
