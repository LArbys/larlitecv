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


unihackalgo = larlitecv.UnipolarHackAlgo()
applyhack = std.vector("int")(3,0)
applyhack[1] = 1
hackthresh = std.vector("float")(3,-10)
img_v = unihackalgo.hackUnipolarArtifact( img_v, applyhack, hackthresh )

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

segmentalgo = larlitecv.Segment3DAlgo()

# ANODE EXAMPLES
#tick_high = 3402
#tick_low  = 3342
#tick_high = 4166
#tick_low  = 4086
# CATHODE EXAMPLES
cathode_adjust = 48.0
#tick_high = 7550
#tick_high = 6332
tick_high = 4566
tick_high -= cathode_adjust
tick_low  = tick_high - 60.0
threshold = 10.0
hitwidth  = 1
frac_good = 0.95

row_high = meta.row( tick_low )
row_low  = meta.row( tick_high  )

hit_high_v = std.vector( "vector<int>" )()
hit_low_v  = std.vector( "vector<int>" )()

for p in range(0,3):
  planehits_hi = segmentalgo.findHits( img_v.at(p), row_high, threshold, hitwidth )
  planehits_lo = segmentalgo.findHits( img_v.at(p), row_low,  threshold, hitwidth )  
  #print "Plane ",p," hits: ",planehits.size()," [",
  #for hid in range(planehits.size()):
  #  print planehits.at(hid)," ",
  #print "]"
  hit_high_v.push_back( planehits_hi )
  hit_low_v.push_back( planehits_lo )

# segment 2d
seg2d_v = std.vector( "vector<larlitecv::Segment2D_t>" )()
for p in range(0,3):
  plane_seg2d_v = segmentalgo.make2DSegments( img_v.at(p), badch_v.at(p), row_low, hit_low_v.at(p), row_high, hit_high_v.at(p), threshold, 1, frac_good )
  print "plane ",p," segments: ",plane_seg2d_v.size(), " hits-low: ",hit_low_v.at(p).size()," hits-high: ",hit_high_v.at(p).size()
  seg2d_v.push_back( plane_seg2d_v )

seg3d_v = std.vector("larlitecv::Segment3D_t")()
threshold_v = std.vector("float")(3,threshold)
segmentalgo.combine2Dinto3D( seg2d_v, img_v, badch_v, hitwidth, threshold_v, frac_good, seg3d_v )

print "Number of Segment-3D: ",seg3d_v.size()

seg3d_all_v = segmentalgo.find3DSegments( img_v, badch_v, row_low, row_high, threshold_v, hitwidth )

print "Number of Segment-3D: ",seg3d_all_v.size()

  

# output image
for p in range(0,3):
  imgarr = larcv.as_ndarray( img_v.at(p) )
  imgarr = np.transpose( imgarr )
  imgarr[ imgarr<threshold ] = 0
  imgarr[ imgarr>255+threshold ] = 255+threshold
  imgarr -= threshold
  cvimg  = np.zeros( (imgarr.shape[0], imgarr.shape[1], 3), dtype=np.int32 )
  cvimg[:,:,0] = imgarr
  cvimg[:,:,1] = imgarr
  cvimg[:,:,2] = imgarr

  seg_v = seg2d_v.at(p)
  for sid in range(seg_v.size()):
    seg = seg_v.at(sid)
    cv2.line( cvimg, (seg.col_low,seg.row_low), (seg.col_high,seg.row_high), (0,255,0,125), 1 )

  # mark 3d segments projected into plane
  for sid in range(seg3d_all_v.size()):
    seg3d = seg3d_all_v.at(sid)
    pixels_v = std.vector("larcv::Pixel2DCluster")()
    s = std.vector("float")(3)
    e = std.vector("float")(3)
    for i in range(3):
      s[i] = seg3d.start[i]
      e[i] = seg3d.end[i]
    thresholds = std.vector("float")(3,10.0)
    neighborhoods = std.vector("int")(3,1)
    larcv.UBWireTool.pixelsAlongLineSegment( s, e, img_v, thresholds, neighborhoods, 0.3, pixels_v )
    for pid in range( pixels_v.at(p).size() ):
      pix = pixels_v.at(p).at(pid)
      x1 = (pix.X(),pix.Y())
      x2 = (x1[0]+1,x1[1]-1)
      cv2.rectangle(cvimg, x1,x2,(255,0,255,20),1)

  # mark hits
  hits_hi = hit_high_v.at(p)
  hits_lo = hit_low_v.at(p)  
  print p,imgarr.shape
  for hid in range(hits_hi.size()):
    hitcenter = hits_hi.at(hid)
    x1 = (hitcenter-hitwidth,row_high+1)
    x2 = (hitcenter+hitwidth,row_high-1)
    cv2.rectangle( cvimg, x1, x2, (0,0,255,125), 1 )
  for hid in range(hits_lo.size()):
    hitcenter = hits_lo.at(hid)
    x1 = (hitcenter-hitwidth,row_low+1)
    x2 = (hitcenter+hitwidth,row_low-1)
    cv2.rectangle( cvimg, x1, x2, (255,0,0,125), 1 )
      
    
  cv2.imwrite( "segimg2d_p%d.png"%(p), cvimg )
  
print "[ENTER] to finish."
raw_input()
