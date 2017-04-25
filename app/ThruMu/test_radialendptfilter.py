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

# Get Endpts
# Load some end points\n,
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

# make cvimage
# output image
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

for c in range(meta.cols()):
  for p in range(3):
    if badch_v[p].pixel(1,c)>0:
      cvimg[:,c,p] = 20
      #cvimg[:,c,1] = 5
      #cvimg[:,c,2] = 5


# SETUP THE ALGO
radialfilter = larlitecv.RadialEndpointFilter()
cfg = larlitecv.RadialEndpointFilterConfig()
cfg.segment_radius = 8.0
cfg.segment_frac_w_charge = 0.5


# RUN OVER LOOPS      
for boundary in ["topspacepts","botspacepts","anodepts","cathodepts"]:
  npts = ev_boundary_pts[boundary].Pixel2DArray(0).size()

  #if boundary not in ["topspacepts"]:
  #  continue
  
  for endptid in range(npts):
    #if endptid not in [6]:
    #  continue

    # run 3D
    row = int(ev_boundary_pts[boundary].Pixel2DArray(0).at(endptid).Y())
    tick = meta.pos_y(row)
    cols = std.vector("int")(3)
    poszy = std.vector("float")()
    for p in range(0,3):
      cols[p] = int(meta.pos_x( ev_boundary_pts[boundary].Pixel2DArray(p).at(endptid).X() ))
      pos3d = std.vector("float")(3)      
      triarea = rt.Double(0)
      crosses = rt.Long(0)
    larcv.UBWireTool.wireIntersection( cols, poszy, triarea, crosses )
    pos3d[0] = ( tick-3200.0 )*(larutil.LArProperties.GetME().DriftVelocity()*0.5)
    pos3d[1] = poszy[1];
    pos3d[2] = poszy[0];
    print "[ ",boundary, " #",endptid," pos=(",pos3d[0],",",pos3d[1],",",pos3d[2],")"," ]"

    num_seg3d = rt.Long(0)
    bad = False
    #try:
    isonline = radialfilter.isWithinStraightSegment( pos3d, img_v, badch_v, cfg, num_seg3d )
    #except:
    #  isonline = False
    #  bad = True
          
    # output image
    if not bad:
      if num_seg3d==0:
        color = (255,255,255)
      else:
        if not isonline:
          color = (255,125,0)
        elif isonline:
          color = (255,0,255)
    else:
      color = (0,0,255)
          

    if boundary in ["topspacepts","botspacepts"]:
      for p in range(0,3):
        cv2.circle( cvimg, (cols[p],row), 5, color, 1 )        
        cv2.putText( cvimg, "%s%dp%d"%(boundary[0],endptid,p), (cols[p]+5, row+3), cv2.FONT_HERSHEY_SIMPLEX, 0.3, color)
    else:
      for p in range(0,3):
        cv2.rectangle( cvimg, (cols[p]-2,row-2), (cols[p]+2,row+2), color, 1 )        
        cv2.putText( cvimg, "%s%dp%d"%(boundary[0],endptid,p), (cols[p]+5, row+3), cv2.FONT_HERSHEY_SIMPLEX, 0.3, color)
      

          
cv2.imwrite( "radialfilter_test.png", cvimg )

    
  
print "[ENTER] to finish."
raw_input()
