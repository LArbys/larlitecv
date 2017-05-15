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

testfile = "../ChargeSegmentAlgos/test/tagger_anaout_larcv.root"
img_producer = "modimg"
badch_producer = "gapchs"
makebadch = False
run_2d = False
label_badchs = True

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
    print "Number of ",endptname,": ",ev_boundary_pts[endptname].Pixel2DArray(0).size()
    for ipt in range(ev_boundary_pts[endptname].Pixel2DArray(0).size()):
        print endptname," #%d"%(ipt),": ",
        print "tick=",meta.pos_y( ev_boundary_pts[endptname].Pixel2DArray(0).at(ipt).Y())," ",
        for p in range(0,3):
            print ev_boundary_pts[endptname].Pixel2DArray(p).at(ipt).X()," ",
        print


# pause
print "Go on"
raw_input()

# make cvimage for drawing
threshold = 10.0
maxrgb = 255
imgshp = ( meta.rows(), meta.cols(), 3 )
cvimg  = np.zeros( imgshp )
for p in range(0,3):
  imgarr = larcv.as_ndarray( img_v.at(p) )
  imgarr = np.transpose( imgarr )
  imgarr[ imgarr<threshold ] = 0
  imgarr[ imgarr>maxrgb+threshold ] = maxrgb+threshold
  imgarr -= threshold

  if p==0:
    cvimg[:,:,2] += imgarr
  elif p==1:
    cvimg[:,:,1] += imgarr
  elif p==2:
    cvimg[:,:,0] += imgarr
    cvimg[:,:,1] += imgarr

  cvimg[ cvimg>maxrgb ] = 255

# badchannels
if label_badchs:
  for c in range(meta.cols()):
    for p in range(3):
      if badch_v[p].pixel(1,c)>0:
        cvimg[:,c,p] += 20
        #cvimg[:,c,1] = 5
        #cvimg[:,c,2] = 5
  cvimg[ cvimg>maxrgb ] = 255


# SETUP THE ALGO
cfg = larlitecv.ThruMuFoxTrotConfig()
cfg.use_thrumu_lead = False
cfg.foxtrotalgo_cfg.segment_radius = 8.0
cfg.foxtrotalgo_cfg.segment_frac_w_charge = 0.5
cfg.foxtrotalgo_cfg.num_step_attempts = 4
cfg.foxtrotalgo_cfg.hit_neighborhood = 1
cfg.foxtrotalgo_cfg.verbosity = 3

foxalgo = larlitecv.ThruMuFoxTrot( cfg )

typedict = { "topspacepts": larlitecv.kTop,
             "botspacepts": larlitecv.kBottom,
             "anodepts":larlitecv.kAnode,
             "cathodepts":larlitecv.kCathode,
             "imgendpts":larlitecv.kImageEnd }

# RUN OVER LOOPS

# EASY TRACK
#start_pt_type = "topspacepts"
#start_pt_idx = 3
#end_pt_type = "cathodepts"
#end_pt_idx = 2

# LOT'SO GAPS. Also, scattering point somewhere in the middle
#start_pt_type = "anodepts"
#start_pt_idx = 3
#end_pt_type = "cathodepts"
#end_pt_idx = 0

# Has small gap
#start_pt_type = "botspacepts"
#start_pt_idx = 4
#end_pt_type = "cathodepts"
#end_pt_idx = 4

# Horizontal Gap + dead channels
start_pt_type = "topspacepts"
start_pt_idx  = 4
end_pt_type   = "botspacepts"
end_pt_idx    = 5

start_pix_v = std.vector('larcv::Pixel2D')()
end_pix_v = std.vector('larcv::Pixel2D')()
for p in range(3):
  start_pix_v.push_back( ev_boundary_pts[start_pt_type].Pixel2DArray(p).at(start_pt_idx) )
  end_pix_v.push_back(   ev_boundary_pts[end_pt_type].Pixel2DArray(p).at(end_pt_idx) )

start_sp = larlitecv.Pixel2SpacePoint( start_pix_v, typedict[start_pt_type], img_v.front().meta() )
end_sp   = larlitecv.Pixel2SpacePoint( end_pix_v,   typedict[end_pt_type], img_v.front().meta() )

# reconstitute boundary space point

# run 3D
track = foxalgo.findThruMuFoxTrack( start_sp, end_sp, img_v, badch_v )

for istep in range(track.size()):
  step =  track[istep]
  cols = []
  pos = std.vector("double")(3,0)
  for i in range(3):
    pos[i] = step.pos()[i]
  tick = step.pos()[0]/(larutil.LArProperties.GetME().DriftVelocity()*0.5)+3200.0
  row  = img_v.front().meta().row( tick )
  if foxalgo.getStepTypes()[istep]==0:
    color=(255,0,255)
  else:
    color=(0,255,255)
  for p in range(3):
    wire = larutil.Geometry.GetME().WireCoordinate( pos, p )
    col = img_v.front().meta().col( wire )
    cv2.circle( cvimg, (col,row), 1, color, 1 )

cv2.imwrite( "thrumufoxtrot_test.png", cvimg )



print "[ENTER] to finish."
raw_input()
