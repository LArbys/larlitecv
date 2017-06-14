import os,sys
import ROOT as rt
from ROOT import std
from larcv import larcv
from larlitecv import larlitecv
from math import fabs
import numpy as np

try:
  import cv2
  has_cv = True
  print "Has OpenCV"
except:
  has_cv = False
  print "No OpenCV"  

rt.gStyle.SetOptStat(0)

ioman = larcv.IOManager( larcv.IOManager.kREAD );
ioman.add_in_file( "/mnt/disk1/production/mcc8_test/bnb_nue_intrinsic/larcv_0000_0099.root" );
ioman.initialize()
ioman.read_entry(0)
ev_img      = ioman.get_data( larcv.kProductImage2D,  "tpc" )
ev_chstatus = ioman.get_data( larcv.kProductChStatus, "tpc" )
img_v   = ev_img.Image2DArray()

# get meta
meta = ev_img.Image2DArray().front().meta()

# make gapchs
emptyalgo = larlitecv.EmptyChannelAlgo()
badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, ev_chstatus )
gapchimgs_v = emptyalgo.findMissingBadChs( img_v, badch_v, 5, 200 );
# combine with badchs
for p in range(0,3):
  gapchimg = gapchimgs_v.at(p);
  gapchimg += badch_v.at(p);
badch_v = gapchimgs_v

print "Ready."

# ==============================================
# ENTRY 0
# example where one track needs to get fixed

runhack = std.vector("int")(3,0)
runhack[1] = 1
thresh = std.vector("float")(3,-10)

unihackalgo    = larlitecv.UnipolarHackAlgo()
#for c in range(img_v.at(1).meta().cols()):
#for c in range(1604,1605):
#  unihack_pulses = std.vector("larlitecv::UnipolarROI_t")()  
#  unihackalgo.scanForPulse( c, img_v.at(1), -5, unihack_pulses )
#  print "number of Unipolar ROIs. in col=%d:"%(c),unihack_pulses.size()

hacked_v = unihackalgo.hackUnipolarArtifact( img_v, runhack, thresh )

for p in range(hacked_v.size()):
  hacked = hacked_v.at(p)
  imgarr = larcv.as_ndarray( hacked_v.at(p) )
  imgarr = np.transpose( imgarr )
  imgarr[ imgarr<10 ] = 0
  imgarr[ imgarr>255 ] = 255
  cvimg = np.zeros( (imgarr.shape[0],imgarr.shape[1], 3 ), dtype=np.int32 )
  cvimg[:,:,0] = imgarr
  cvimg[:,:,1] = imgarr
  cvimg[:,:,2] = imgarr
  cv2.imwrite( "hackedimg2d_p%d.png"%(p), cvimg )
  
print "[ENTER] to finish."
raw_input()
