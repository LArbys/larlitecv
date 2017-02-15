import os,sys
import ROOT as rt
from ROOT import std
from larcv import larcv
from larlitecv import larlitecv
from math import fabs


# we test bezier 2D fit
# because we will have parameterized curve, we can in principle piece together a 3D curve after combining 3 (or 2) plane information

rt.gStyle.SetOptStat(0)
padding = 20
badchval = -20

algoconfig = larlitecv.Linear3DFitterConfig()
algoconfig.neighborhood_posttick = 10

# ------------------------------------------------------------------------------
# SETUP TEST PATH

# example of track that has to bend and not jump using bad channels
#entry = 0
#start_tick  = 3348
#end_tick    = 2412
#start_wires = [202,477,11]
#end_wires   = [696,34,60]

# straight path that ends in dead region
#entry = 2
#start_tick = 4884
#start_wires = [465,1122,913]
#end_tick = 6282
#end_wires = [851,247,426]

# example of spurious connection
#entry = 2
#start_tick = 6762
#start_wires = [1061,414,802]
#end_tick = 6360 # bad ending
#end_wires = [585,1092,1005] # bad ending
#end_tick = 6054 # correct ending
#end_wires = [565,1220,1114] # correct ending

#entry = 2
#start_tick = 4792
#start_wires = [1392,2006,2726]
#end_tick = 3000
#end_wires = [1927,1324,2580]

# tuning end-point rejection efficiency
# example of determining what should be connect
#entry = 2
#start_tick = 7176
#start_wires = [339,1008,675]
#end_tick = 8046 # correct end#
#end_wires = [866,206,398] # correct end
#end_tick  = 7230  # wrong end point in middle of track
#end_wires = [374,960,662] # wrong end point
#end_tick = []

entry = 2
start_tick = 6054
start_wires = [565,1220,1114]
end_tick    = 6762
end_wires   = [1061,414,802]

#entry = 0
#start_tick = 5046
#start_wires = [1272,1670,2269]
#end_tick = 7986
#end_wires = [1689,1085,2101]

#entry = 1
#start_tick = 7890
#start_wires = [1001,1597,1928]
#end_tick = 6084
#end_wires = [1618,1015,1960]

#entry = 5
#start_tick = 4176
#start_wires = [2311,1806,3445]
#end_tick = 2418
#end_wires = [1630,2182,3135]

#entry = 9
#start_tick = 4680 
#start_wires = [1081,1742,2152]
#end_tick = 4890
#end_wires = [2176,1506,3013]

# ------------------------------------------------------------------------------
# first let's prepare data

ioman = larcv.IOManager( larcv.IOManager.kREAD );
ioman.add_in_file( "bin/output_larcv_testextbnb_run00.root" );
ioman.initialize()
ioman.read_entry(entry)
ev_img   = ioman.get_data( larcv.kProductImage2D, "modimgs" )
ev_badchs = ioman.get_data( larcv.kProductImage2D, "gapchs" )
img_v   = ev_img.Image2DArray()
badch_v = ev_badchs.Image2DArray()

# get meta
meta = ev_img.Image2DArray().front().meta()

# make gapchs
emptyalgo = larlitecv.EmptyChannelAlgo()
gapchimgs_v = emptyalgo.findMissingBadChs( img_v, badch_v, 5, 200 );
# combine with badchs
for p in range(0,3):
  gapchimg = gapchimgs_v.at(p);
  gapchimg += badch_v.at(p);
badch_v = gapchimgs_v

# ------------------------------------------------------------------------------
# fill 2D histograms with orginal image containing track
himgs = []
for p in range(0,3):
	low_wire = start_wires[p]
	high_wire = end_wires[p]
	if low_wire>high_wire:
		low_wire = end_wires[p]
		high_wire = start_wires[p]
	high_wire += padding*meta.pixel_width()
	low_wire  -= padding*meta.pixel_width()
	if low_wire<0:
		low_wire = 0
	print "wire range: ",low_wire,high_wire

	low_tick = start_tick
	high_tick = end_tick
	if low_tick>high_tick:
		low_tick = end_tick
		high_tick = start_tick
	low_tick -= padding*meta.pixel_height()
	if low_tick<=meta.min_y():
		low_tick = meta.min_y()+1
	high_tick += padding*meta.pixel_height()
	if high_tick>=meta.max_y():
		high_tick = meta.max_y()-1

	nxbins = int( fabs(high_wire-low_wire)/meta.pixel_width() )
	nybins = int( (high_tick-low_tick)/meta.pixel_height() )

	print "bins x=",nxbins," y=",nybins
	print "endy=",low_tick + meta.pixel_height()*nybins
	print "endx=",low_wire + meta.pixel_width()*nxbins

	himg = rt.TH2D("himg_p%d"%(p),"himg", nxbins, low_wire, low_wire+meta.pixel_width()*nxbins, nybins, low_tick, low_tick + meta.pixel_height()*nybins )

	badchs = badch_v.at(p)	
	for x in range(0,nxbins):
		for y in range(0,nybins):
			w = himg.GetXaxis().GetBinLowEdge(x+1)
			t = himg.GetYaxis().GetBinLowEdge(y+1)
			if t<= meta.min_y():
				t = meta.min_y()+1;
			elif t>=meta.max_y():
				t = meta.max_y()-1;
			if w>=meta.cols():
				w = meta.cols()-1
			col = meta.col( w )
			row = meta.row( t )
			val = img_v.at(p).pixel( row, col )
			if badchs.pixel(row,col)>0:
				himg.SetBinContent(x+1,y+1,badchval)
				continue			
			if val<10:
				continue
			himg.SetBinContent( x+1, y+1, val )
	himgs.append(himg)

# Algo setup

algo = larlitecv.Linear3DFitter(algoconfig)
start_row = img_v.front().meta().row( start_tick )
end_row   = img_v.front().meta().row( end_tick )
start_cols = std.vector("int")(3,0)
end_cols   = std.vector("int")(3,0)
for i in range(0,3):
	start_cols[i] = img_v.front().meta().col( start_wires[i] )
	end_cols[i]   = img_v.front().meta().col( end_wires[i] )	
pointlist = algo.findpath( img_v, badch_v, start_row, end_row, start_cols, end_cols )
print "Number of points in list: ",pointlist.size()
print "Fraction w/ charge: ",pointlist.fractionHasChargeWith3Planes()
print "Fraction w/ badch: ",pointlist.fractionHasBadChOn3Planes()
print "Fraction w/ no charge: ",pointlist.fractionHasNoChargeOn3Planes()
print "Fraction w/ majority charge: ",pointlist.fractionHasChargeOnMajorityOfPlanes()
print "Fraction good: ",pointlist.fractionGood()

start_ext = larlitecv.PointInfoList()
end_ext   = larlitecv.PointInfoList()
algo.getTrackExtension( pointlist, img_v, badch_v, 10.0, start_ext, end_ext )
print "start extension: ",start_ext.size()
print "  ngood: ",start_ext.num_pts_good
print "  nallcharge: ",start_ext.num_pts_w_allcharge
print "  nmajcharge: ",start_ext.num_pts_w_majcharge
print "  nallbad: ",start_ext.num_pts_w_allbadch
print "end extension: ",end_ext.size()
print "  ngood: ",end_ext.num_pts_good
print "  nallcharge: ",end_ext.num_pts_w_allcharge
print "  nmajcharge: ",end_ext.num_pts_w_majcharge
print "  nallbad: ",end_ext.num_pts_w_allbadch


# ----------------------------------------------------------------
# Visualize

npts = pointlist.size()
gpaths = []
gbadch = []
gwq    = []
for p in range(0,3):
	g = rt.TGraph(npts)
	g.SetMarkerStyle(20)
	for pt in range(0,npts):
		g.SetPoint(pt, pointlist.at(pt).wire_id.at(p), pointlist.at(pt).tick  )
	g.SetMarkerColor( rt.kBlack )	
	gpaths.append(g)

	g = rt.TGraph( pointlist.num_pts_good )
	ipt = 0
	for pt in range(0,npts):
		if pointlist.at(pt).goodpoint==True:
			g.SetPoint(ipt,pointlist.at(pt).wire_id.at(p),pointlist.at(pt).tick)
			ipt += 1
	g.SetMarkerColor( rt.kRed )
	g.SetMarkerStyle(20)
	gbadch.append(g)

	g = rt.TGraph( pointlist.num_pts_w_allcharge )
	ipt = 0
	for pt in range(0,npts):
		if pointlist.at(pt).planeswithcharge==3 and pointlist.at(pt).goodpoint==True:
			g.SetPoint(ipt,pointlist.at(pt).wire_id.at(p),pointlist.at(pt).tick)
			ipt += 1
	g.SetMarkerColor( rt.kYellow )
	g.SetMarkerStyle(20)
	gwq.append(g)

c = rt.TCanvas("c1","c1",1200,400)
c.Divide(3,1)

for i in range(0,3):
	c.cd(i+1)
	himgs[i].Draw("COLZ")
	gpaths[i].Draw("LP")
	gbadch[i].Draw("P")
	gwq[i].Draw("P")

c.Update()
print "FINISHED. [Enter] to quit."
raw_input()



