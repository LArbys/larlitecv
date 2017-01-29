import os,sys
import ROOT as rt
from ROOT import std
from larcv import larcv
from larlitecv import larlitecv
from math import fabs

rt.gStyle.SetOptStat(0)

ioman = larcv.IOManager( larcv.IOManager.kREAD );
ioman.add_in_file( "bin/output_larcv_testextbnb_run00.root" );
ioman.initialize()
ioman.read_entry(0)
ev_img   = ioman.get_data( larcv.kProductImage2D, "modimgs" )
ev_badchs = ioman.get_data( larcv.kProductImage2D, "gapchs" )
img_v   = ev_img.Image2DArray()
badch_v = ev_badchs.Image2DArray()

# get meta
meta = ev_img.Image2DArray().front().meta()

# example where one track needs to get fixed
#start_tick = 5634
#end_tick   = 3876
#start_wires = [731,1398,1456]
#end_wires   = [1795,1124,2246]

#start_tick = 3978
#end_tick   = 2418
#start_wires = [633,1295,1257]
#end_wires   = [1171,838,1337]

# track that has to bend
start_tick  = 3348
end_tick    = 2412
start_wires = [202,477,11]
end_wires   = [696,34,60]

# example where end point is in dead region
use_bad_chs = False

# setup algo
config  = larlitecv.AStar3DAlgoConfig()
config.astar_threshold.resize(3,50.0)
config.astar_threshold[2] = 200.0
config.astar_neighborhood.resize(3,2)
config.astar_start_padding = 3
config.astar_end_padding = 3
config.lattice_padding = 10

algo = larlitecv.AStar3DAlgo( config )
algo.setVerbose(2)

# attempt compress until under 100 per dimension
img_compress_v   = std.vector("larcv::Image2D")()
badch_compress_v = std.vector("larcv::Image2D")()
for p in range(0,3):
	img_cpy   = larcv.Image2D( img_v.at(p) )
	badch_cpy = larcv.Image2D( badch_v.at(p) )
	img_cpy.compress( img_v.at(p).meta().rows()/4, img_v.at(p).meta().cols()/4 )
	badch_cpy.compress( img_v.at(p).meta().rows()/4, img_v.at(p).meta().cols()/4 )	
	img_compress_v.push_back( img_cpy )
	badch_compress_v.push_back( badch_cpy )

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

path = algo.findpath( img_compress_v, badch_compress_v, start_row, end_row, start_cols, end_cols )

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

	nxbins = int( fabs(high_wire-low_wire)/meta.pixel_width() )
	nybins = int( (high_tick-low_tick)/meta.pixel_height() )

	print "bins x=",nxbins," y=",nybins
	print "endy=",low_tick + meta.pixel_height()*nybins
	print "endx=",low_wire + meta.pixel_width()*nxbins

	himg = rt.TH2D("himg_p%d"%(p),"himg", nxbins, low_wire, low_wire+meta.pixel_width()*nxbins, nybins, low_tick, low_tick + meta.pixel_height()*nybins )

	badchs = badch_v.at(p)	
	for x in range(0,nxbins):
		for y in range(0,nybins):
			col = meta.col( himg.GetXaxis().GetBinLowEdge(x+1) )
			row = meta.row( himg.GetYaxis().GetBinLowEdge(y+1) )
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

for i in range(0,path.size()):
	node = path.at(i)
	for p in range(0,3):
		x = meta_compress.pos_x( node.cols[p] )
		y = meta_compress.pos_y( node.row )
		gpaths[p].SetPoint( i, x, y )

c = rt.TCanvas("c","c",1400,1200)
c.Divide(3,3)
c.Draw()

for p in range(0,3):
	c.cd( 3*p+1)
	himgs[p].Draw("COLZ")
	gpaths[p].Draw("LP")
	c.cd(3*p+2)
	#hfscores[p].Draw("COLZ")
	c.cd(3*p+3)
	#hgscores[p].Draw("COLZ")

c.Update()

raw_input()
