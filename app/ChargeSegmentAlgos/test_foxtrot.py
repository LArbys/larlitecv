import os,sys
import ROOT as rt
from ROOT import std
from larcv import larcv
from larlitecv import larlitecv
from ROOT import larutil
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

testfile = "../ChargeSegmentAlgos/test/tagger_anaout_larcv_seg.root"
img_producer = "modimg"
badch_producer = "gapchs"
makebadch = False
run_2d = False

ioman = larcv.IOManager( larcv.IOManager.kREAD );
ioman.add_in_file( testfile )
ioman.initialize()
ioman.read_entry(0)
ev_img      = ioman.get_data( larcv.kProductImage2D,  img_producer )
img_v   = ev_img.Image2DArray()

unihackalgo = larlitecv.UnipolarHackAlgo()
applyhack = std.vector("int")(3,0)
applyhack[1] = 1
hackthresh = std.vector("float")(3,-10)
img_v = unihackalgo.hackUnipolarArtifact( img_v, applyhack, hackthresh )

# get meta
meta = ev_img.Image2DArray().front().meta()

# make gapchs
if makebadch:
  ev_chstatus = ioman.get_data( larcv.kProductChStatus, "tpc" )
  emptyalgo = larlitecv.EmptyChannelAlgo()
  badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, ev_chstatus )
  gapchimgs_v = emptyalgo.findMissingBadChs( img_v, badch_v, 5, 200 );
  # combine with badchs
  for p in range(0,3):
    gapchimg = gapchimgs_v.at(p);
    gapchimg += badch_v.at(p);
    badch_v = gapchimgs_v
else:
  ev_badch_v = ioman.get_data( larcv.kProductImage2D, badch_producer )
  badch_v = ev_badch_v.Image2DArray()


print "Images Ready."

# Display the end points,
endptnames = ["topspacepts","botspacepts","upspacepts","downspacepts","anodepts","cathodepts","imgendpts"]
ev_boundary_pts = {}
for endptname in endptnames:
    ev_boundary_pts[endptname] = ioman.get_data(larcv.kProductPixel2D, endptname )
    for ipt in range(ev_boundary_pts[endptname].Pixel2DArray(0).size()):
        print endptname," #%d"%(ipt),": ",
        print "tick=",meta.pos_y( ev_boundary_pts[endptname].Pixel2DArray(0).at(ipt).Y())," ",
        for p in range(0,3):
            print ev_boundary_pts[endptname].Pixel2DArray(p).at(ipt).X()," ",
        print

# make cvimage for drawing
threshold = 10.0
imgshp = ( meta.rows(), meta.cols(), 3 )
cvimg  = np.zeros( imgshp )
for p in range(0,3):
  imgarr = larcv.as_ndarray( img_v.at(p) )
  imgarr = np.transpose( imgarr )
  imgarr[ imgarr<threshold ] = 0
  imgarr[ imgarr>255+threshold ] = 255+threshold
  imgarr -= threshold

  if p==0:
    cvimg[:,:,2] += imgarr
  elif p==1:
    cvimg[:,:,1] += imgarr
  elif p==2:
    cvimg[:,:,0] += imgarr
    cvimg[:,:,1] += imgarr

  cvimg[ cvimg>255 ] = 255

# badchannels
for c in range(meta.cols()):
  for p in range(3):
    if badch_v[p].pixel(1,c)>0:
      cvimg[:,c,p] = 20
      #cvimg[:,c,1] = 5
      #cvimg[:,c,2] = 5


# SETUP THE ALGO
cfg = larlitecv.FoxTrotTrackerAlgoConfig();
cfg.segment_radius = 8.0
cfg.segment_frac_w_charge = 0.5
cfg.num_step_attempts = 4
cfg.verbosity = 3

foxalgo = larlitecv.FoxTrotTrackerAlgo( cfg )

typedict = { "topspacepts": larlitecv.kTop,
             "botspacepts": larlitecv.kBottom,
             "anodepts":larlitecv.kAnode,
             "cathodepts":larlitecv.kCathode }

# RUN OVER LOOPS
for boundary in ["topspacepts","botspacepts","anodepts","cathodepts"]:
  npts = ev_boundary_pts[boundary].Pixel2DArray(0).size()

  if boundary not in ["botspacepts"]:
     continue

  for endptid in range(npts):
    if endptid not in [4]:
      continue

    # reconstitute boundary space point
    pix_v = std.vector('larcv::Pixel2D')()
    for p in range(3):
      pix_v.push_back( ev_boundary_pts[boundary].Pixel2DArray(p).at(endptid) )
    sp = larlitecv.Pixel2SpacePoint( pix_v, typedict[boundary], img_v.front().meta() )

    # run 3D
    track = foxalgo.followTrack( img_v, badch_v, badch_v, sp )

    for istep in range(track.size()):
      step =  track[istep]
      cols = []
      pos = std.vector("double")(3,0)
      for i in range(3):
        pos[i] = step.pos()[i]
      tick = step.pos()[0]/(larutil.LArProperties.GetME().DriftVelocity()*0.5)+3200.0
      row  = img_v.front().meta().row( tick )
      for p in range(3):
        wire = larutil.Geometry.GetME().WireCoordinate( pos, p )
        col = img_v.front().meta().col( wire )
        cv2.circle( cvimg, (col,row), 1, (255,0,255), 1 )


cv2.imwrite( "foxtrot_test.png", cvimg )



print "[ENTER] to finish."
raw_input()
