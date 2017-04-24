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

testfile = "test/tagger_anaout_larcv_seg.root"
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


# SETUP THE ALGO
radialalgo = larlitecv.RadialSegmentSearch()

# ==============================================
# ENTRY 0
# example where one track needs to get fixed

# Test point
"#topspacepts  #6 :  tick= 6120.0   214   827   369  \n",
#boundary_type = "botspacepts"
#endptid = 3

def segment_box( boundary_type, endptid, radius ):
  row = int(ev_boundary_pts[boundary_type].Pixel2DArray(0).at(endptid).Y())
  tick = meta.pos_y(row)
  cols = std.vector("int")(3)
  poszy = std.vector("float")()
  for p in range(0,3):
    cols[p] = int(meta.pos_x( ev_boundary_pts[boundary_type].Pixel2DArray(p).at(endptid).X() ))
  pos3d = std.vector("float")(3)      
  triarea = rt.Double(0)
  crosses = rt.Long(0)
  larcv.UBWireTool.wireIntersection( cols, poszy, triarea, crosses )
  pos3d[0] = ( tick-3200.0 )*(larutil.LArProperties.GetME().DriftVelocity()*0.5)
  pos3d[1] = poszy[1];
  pos3d[2] = poszy[0];
  print boundary_type, " #",endptid," pos=(",pos3d[0],",",pos3d[1],",",pos3d[2],")",
  
  plane_points = {}
  for p in range(0,3):
    radial_points = radialalgo.findIntersectingChargeClusters( img_v[p], badch_v[p], pos3d, radius, 10.0 )
    print "  plane ",p," radial points: ",radial_points.size()
    plane_points[p] = radial_points
    for iradpix in range(radial_points.size()):
      radhit = plane_points[p].at(iradpix)
      if radhit.max_idx<0 or radhit.max_idx>=radial_points.size():
        continue
      radpix = radhit.pixlist[ radhit.max_idx ]
      print " tick=", meta.pos_y(radpix.Y())," max=",radpix.X(),
      for ipix in range( radhit.pixlist.size() ):
        radpix = radhit.pixlist[ ipix ]
        print "(%d,%d)"%(meta.pos_y(radpix.Y()),radpix.X())," ",
      print

  # output image
  for p in range(0,3):
    cv2.circle( cvimg, (cols[p],row), 5, (255,255,255), 1 )

    for iradpix in range(plane_points[p].size()):
      radhit = plane_points[p].at(iradpix)
      radpix = radhit.pixlist[ radhit.max_idx ]
      cv2.circle( cvimg, (radpix.X(),radpix.Y()), 2, (255,0,0), -1 )
  return plane_points, pos3d


def make2dseg( plane_points ):
  return
      
for boundary in ["topspacepts","botspacepts"]:
  npts = ev_boundary_pts[boundary].Pixel2DArray(0).size()

  if boundary not in ["topspacepts"]:
    continue
  
  for endptid in range(npts):
    if endptid not in [6]:
      continue

    if run_2d:
      plane_radhits, pos3d = segment_box( boundary, endptid, 10.0 )
      plane_seg2d = {}
      for p in range(3):
        plane_seg2d[p] = radialalgo.make2Dsegments( img_v[p], badch_v[p], plane_radhits[p], pos3d, 10.0, 1, 0.8 )
        print "number of 2d segments: ",plane_seg2d[p].size()
        for iseg in range(plane_seg2d[p].size()):
          seg2d = plane_seg2d[p][iseg]
          cv2.line( cvimg, (seg2d.col_low,seg2d.row_low), (seg2d.col_high,seg2d.row_high), (0,125,255), 1 )
    else:
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
      print boundary, " #",endptid," pos=(",pos3d[0],",",pos3d[1],",",pos3d[2],")"

      pix_thresholds = std.vector("float")(3,10.0)
      seg3d_v = radialalgo.find3Dsegments( img_v, badch_v, pos3d, 10.0, pix_thresholds, 1, 0.8 )
      print "  number of 3d segments=", seg3d_v.size()
          



          
cv2.imwrite( "radial_test.png", cvimg )

    
  
print "[ENTER] to finish."
raw_input()
