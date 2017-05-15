import os,sys
import ROOT as rt
from ROOT import std
from ROOT import Long, Double
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
colors = {0:(255,0,255,255),
          1:(255,255,0,255),
          2:(0,255,255,255)}

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
maxval = 100.0
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

  cvimg[ cvimg>maxval ] = 255


# SETUP THE ALGO
radialalgo = larlitecv.RadialSegmentSearch()

# ==============================================
# ENTRY 0
# example where one track needs to get fixed

# Test points
"#topspacepts  #6 :  tick= 6120.0   214   827   369  \n",
"#topspacepts  #12 :  tick= 3192.0   1257   1911   2497 \n" # horizontal segments

boundary_type = "topspacepts"
endptid = 12

# get 3d point

cm_per_tick = larutil.LArProperties.GetME().DriftVelocity()*0.5
wids = std.vector("int")(3)
for p in range(0,3):
  wids[p] = ev_boundary_pts[boundary_type].Pixel2DArray(p).at(endptid).X()
row = ev_boundary_pts[boundary_type].Pixel2DArray(0).at(endptid).Y()
tick = img_v.front().meta().pos_y( ev_boundary_pts[boundary_type].Pixel2DArray(0).at(endptid).Y() )
x = (img_v.front().meta().pos_y( ev_boundary_pts[boundary_type].Pixel2DArray(0).at(endptid).Y() ) - 3200.0)*cm_per_tick

poszy = std.vector("float")()
crosses = Long(0)
triarea  = Double(0)
larcv.UBWireTool.wireIntersection( wids, poszy, triarea, crosses )
pos3d = std.vector("float")(3)
pos3d[0] = x
pos3d[1] = poszy[1]
pos3d[2] = poszy[0]

print "Test point 3d pos: (",pos3d[0],",",pos3d[1],",",pos3d[2],"), img coords: (",tick,",",wids[0],",",wids[1],",",wids[2],")"
for p in range(0,3):
  cv2.circle( cvimg, (wids[p],row), 1, (255,255,255,255), 1 )
pixel_thresholds = std.vector("float")(3,10.0)

seg3d_v = radialalgo.find3Dsegments( img_v, badch_v, pos3d, 8.0, pixel_thresholds, 1, 1, 0.8, 3 )
print "Number of 3D segments: ",seg3d_v.size()

for p in range(3):
  for ihit in range(radialalgo.m_planehits[p].size()):
    radhit = radialalgo.m_planehits[p].at(ihit)
    cv2.circle( cvimg, (radhit.col(),radhit.row()), 1, colors[p], 1 )


for iseg in range(seg3d_v.size()):
  seg3d = seg3d_v[iseg]
  startpos3d = std.vector("float")(3)
  endpos3d   = std.vector("float")(3)
  for i in range(3):
    startpos3d[i] = seg3d.start[i]
    endpos3d[i] = seg3d.end[i]
  startimgcoord = larcv.UBWireTool.getProjectedImagePixel( startpos3d, img_v.front().meta(), 3 )
  endimgcoord   = larcv.UBWireTool.getProjectedImagePixel( endpos3d, img_v.front().meta(), 3 )

  print "iseg #",iseg,": (",startpos3d[0],",",startpos3d[1],",",startpos3d[2],") -> (",endpos3d[0],",",endpos3d[1],",",endpos3d[2],")"

  for p in range(3):
    cv2.circle( cvimg, (startimgcoord[p+1],startimgcoord[0]), 3, (255,255,255,255), 1 )
    cv2.circle( cvimg, (endimgcoord[p+1],endimgcoord[0]), 3, colors[p], 1 )


cv2.imwrite( "radial_test.png", cvimg )



print "[ENTER] to finish."
raw_input()
