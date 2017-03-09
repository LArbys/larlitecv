import os,sys
import ROOT as rt
from ROOT import std
from larcv import larcv
from larlitecv import larlitecv
from math import fabs

rt.gStyle.SetOptStat(0)

ioman = larcv.IOManager( larcv.IOManager.kREAD );
ioman.add_in_file( "astar_test_file.root" );
#ioman.add_in_file( "bin/output_larcv_testextbnb_run00.root" );
ioman.initialize()
ioman.read_entry(0)
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

# ==============================================
# ENTRY 0
# example where one track needs to get fixed
#start_tick = 5634
#end_tick   = 3876
#start_wires = [731,1398,1456]
#end_wires   = [1795,1124,2246]

#start_tick = 3978
#end_tick   = 2418
#start_wires = [633,1295,1257]
#end_wires   = [1171,838,1337]

# track that has to bend and not jump using bad channels
#start_tick  = 3348
#end_tick    = 2412
#start_wires = [202,477,11]
#end_wires   = [696,34,60]

# example of spinning its wheels for a long time on a false path that passes heuristic
# plane 2 is good path. plane 0 is full of badchs, plane 1 has no path
#[ path-finding for endpoints (15,44) of type (0) -> (5) ]
#  plane=0:  (w,t): (625, 4944) -> (904,5784)
#  plane=1:  (w,t): (1285, 4944) -> (803,5784)
#  plane=2:  (w,t): (1241, 4944) -> (1035,5784)
#  line region test: 0.978417, 0.453222, 1
#start_tick = 4944
#end_tick = 5784
#start_wires = [625,1285,1241]
#end_wires = [904,803,1035]

# example where one can get stuck inside dead region
#start_tick = 7986
#end_tick = 5046
#start_wires = [1689,1085,2101]
#end_wires   = [1272,1670,2269]

# path that crosses a large dead-channel gap 
start_tick  = 4578
start_wires = [239,894,463]
end_tick    = 5730
end_wires   = [1090,490,907]

# snakes through dead channels
#start_tick = 3792
#start_wires = [702,1295,1329]
#end_tick = 3090
#end_wires = [868,1106,1302]

# ==============================================

# ENTRY 1
#start_tick = 7890
#end_tick = 6084
#start_wires = [1001,1597,1928]
#end_wires   = [1618,1015,1960]

# only major track that is missed
#start_tick = 3606
#start_wires = [555,1210,1097]
#end_tick = 5568
#end_wires = [1039,708,1070]

#start_tick = 3540
#start_wires = [514,1053,896]
#end_tick = 5568
#end_wires = [1039,708,1070]

# ==============================================

# ENTRY 2
#start_tick = 2900
#end_tick   = 2418
#start_wires = [1070,1710,2112]
#end_wires  = [1409,1453,2196]

#start_tick = 4050
#start_wires = [1343,2008,2677]
#end_tick = 4644
#end_wires = [1954,1290,2573]

#start_tick = 4884
#start_wires = [465,1122,913]
#end_tick = 6282
#end_wires = [851,247,426]

# example of spurious connection
#start_tick = 6762
#start_wires = [1061,414,802]
#end_tick = 6360 # bad ending
#end_wires = [585,1092,1005] # bad ending
#end_tick = 6054
#end_wires = [565,1220,1114]

#start_tick  = 6054
#start_wires = [565,1220,1114]
#end_tick = 6762
#end_wires = [1061,414,802]


#start_tick = 7176
#start_wires = [339,1008,675]
#end_tick = 8046
#end_wires = [866,206,398]

# example where end point is in dead region
use_bad_chs = False

# setup algo
config  = larlitecv.AStar3DAlgoConfig()
config.accept_badch_nodes = True;
config.astar_threshold.resize(3,50.0)
config.astar_threshold[2] = 50.0
config.astar_neighborhood.resize(3,6)
config.astar_start_padding = 3
config.astar_end_padding = 3
config.lattice_padding = 10
config.min_nplanes_w_hitpixel = 3
config.restrict_path = True
config.path_restriction_radius = 30.0

compress_factor = 4

algo = larlitecv.AStar3DAlgo( config )
algo.setVerbose(2)

# attempt compress until under 100 per dimension
img_compress_v    = std.vector("larcv::Image2D")()
badch_compress_v  = std.vector("larcv::Image2D")()
tagged_v          = std.vector("larcv::Image2D")()
tagged_compress_v = std.vector("larcv::Image2D")()
for p in range(0,3):
	img_cpy   = larcv.Image2D( img_v.at(p) )
	badch_cpy = larcv.Image2D( badch_v.at(p) )
	tag_cpy   = larcv.Image2D( img_v.at(p).meta() )
	tag_cpy.paint(0)
	tagged_v.push_back( tag_cpy )
	img_cpy.compress( img_v.at(p).meta().rows()/compress_factor, img_v.at(p).meta().cols()/compress_factor )
	badch_cpy.compress( img_v.at(p).meta().rows()/compress_factor, img_v.at(p).meta().cols()/compress_factor )	
	tag_cpy.compress( img_v.at(p).meta().rows()/compress_factor, img_v.at(p).meta().cols()/compress_factor )	
	img_compress_v.push_back( img_cpy )
	badch_compress_v.push_back( badch_cpy )
	tagged_compress_v.push_back( tag_cpy )

meta_compress = img_compress_v.front().meta()

start_row  = meta_compress.row( start_tick )
end_row    = meta_compress.row( end_tick )
start_cols = std.vector("int")(3)
end_cols   = std.vector("int")(3)
for p in range(0,3):
	start_cols[p] = meta_compress.col(start_wires[p])
	end_cols[p]   = meta_compress.col(end_wires[p])

# fully replace
img_v = img_compress_v
badch_v = badch_compress_v
meta = meta_compress
goal_reached = rt.Long(0)

path = algo.findpath( img_compress_v, badch_compress_v, tagged_compress_v, start_row, end_row, start_cols, end_cols, goal_reached )

print "pathsize=",path.size()

for p in range(0,path.size()):
	print "node ",p,": ",path.at(p).str()

# fill fscore image
if False:
	hfscores = []
	hgscores = []
	for p in range(0,3):
		fimg = fimgs[p]
		gimg = gimgs[p]
		print "fimg: ",fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y()
		hfscore = rt.TH2D("fscores_p%d"%(p),"fscores", fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y() )
		hgscore = rt.TH2D("gscores_p%d"%(p),"gscores", fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y() )
		for r in range(0,fimg.meta().rows()):
			for c in range(0,fimg.meta().cols()):
				rbin = fimg.meta().row( hfscore.GetYaxis().GetBinLowEdge(r+1)+0.1 )
				cbin = fimg.meta().col( hfscore.GetXaxis().GetBinLowEdge(c+1)+0.1 )
				hfscore.SetBinContent( c+1, r+1, fimg.pixel(rbin,cbin) )
				hgscore.SetBinContent( c+1, r+1, gimg.pixel(rbin,cbin) )		
		hfscores.append(hfscore)
		hgscores.append(hgscore)

# fill original image
himgs = []
for p in range(0,3):
	low_wire = start_wires[p]
	high_wire = end_wires[p]
	if low_wire>high_wire:
		low_wire = end_wires[p]
		high_wire = start_wires[p]
	high_wire += config.lattice_padding*meta.pixel_width()
	low_wire  -= config.lattice_padding*meta.pixel_width()
	if low_wire<0:
		low_wire = 0
	print "wire range: ",low_wire,high_wire

	low_tick = start_tick
	high_tick = end_tick
	if low_tick>high_tick:
		low_tick = end_tick
		high_tick = start_tick
	low_tick -= config.lattice_padding*meta.pixel_height()
	if low_tick<=meta.min_y():
		low_tick = meta.min_y()+1
	high_tick += config.lattice_padding*meta.pixel_height()
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
			col = meta.col( w )
			row = meta.row( t )
			val = img_v.at(p).pixel( row, col )
			if badchs.pixel(row,col)>0:
				himg.SetBinContent(x+1,y+1,-500.0)
				continue			
			if val<10:
				continue
			himg.SetBinContent( x+1, y+1, val )
	himgs.append(himg)


# make path graphs
gpaths = []
for p in range(0,3):
	gpath = rt.TGraph( path.size() )
	gpath.SetMarkerStyle(20)	
	gpaths.append(gpath)

badnodes = {0:[],1:[],2:[]}
for i in range(0,path.size()):
	node = path.at(i)
	for p in range(0,3):
		x = meta_compress.pos_x( node.cols[p] )
		r = node.row
		#if p==2 and r-1>=0:
		#	r += -2
		y = meta_compress.pos_y( r )
		gpaths[p].SetPoint( i, x, y )
		if node.badchnode:
			badnodes[p].append( (x,y) )

gbadnodes = []
for p in range(0,3):
	gpath = rt.TGraph( len(badnodes[p]))
	for n,node in enumerate(badnodes[p]):
		gpath.SetPoint(n,node[0],node[1])
	gpath.SetMarkerStyle(20)
	gpath.SetMarkerColor(rt.kRed)
	gbadnodes.append(gpath)

c = rt.TCanvas("c","c",1400,1200)
c.Divide(3,3)
c.Draw()

for p in range(0,3):
	c.cd( 3*p+1)
	himgs[p].Draw("COLZ")
	gpaths[p].Draw("LP")
	if gbadnodes[p].GetN()>0:
		gbadnodes[p].Draw("P")
	c.cd(3*p+2)
	c.cd(3*p+3)
	#hgscores[p].Draw("COLZ")

c.Update()

raw_input()
