import ROOT as rt
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

# example where one track needs to get fixed
#start_tick = 5634
#end_tick   = 3876
#start_wires = [731,1398,1456]
#end_wires   = [1795,1124,2246]

start_tick = 3978
end_tick   = 2418
start_wires = [633,1295,1257]
end_wires   = [1171,838,1337]

# example where end point is in dead region
use_bad_chs = False

# setup algo
config  = larlitecv.AStarDirAlgoConfig()
config.astar_threshold.resize(3,10.0)
config.astar_neighborhood.resize(3,5)
config.astar_start_padding = 2
config.astar_end_padding = 2
config.image_padding = 20

algo = larlitecv.AStarDirAlgo( config )
algo.setVerbose(0)

# get data
meta = ev_img.Image2DArray().at(2).meta()

paths = []
fimgs = []
gimgs = []

# find paths
for p in range(0,3):
	algo.setBadChImage( ev_badchs.Image2DArray().at(p) )
	path = algo.findpath( ev_img.Image2DArray().at(p), meta.row(start_tick), meta.col(start_wires[p]), meta.row(end_tick), meta.col(end_wires[p]), 10.0, use_bad_chs )
	fimg = algo.getScoreImages().at(0)
	gimg = algo.getScoreImages().at(1)
	print "plane %d path size: %d"%(p,path.size()	)
	paths.append(path)
	fimgs.append(larcv.Image2D(fimg))
	gimgs.append(larcv.Image2D(gimg))


# fill fscore image
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
	high_wire += config.image_padding*meta.pixel_width()
	low_wire  -= config.image_padding*meta.pixel_width()
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

	badchs = ev_badchs.Image2DArray().at(p)	
	for x in range(0,nxbins):
		for y in range(0,nybins):
			col = meta.col( himg.GetXaxis().GetBinLowEdge(x+1) )
			row = meta.row( himg.GetYaxis().GetBinLowEdge(y+1) )
			val = ev_img.Image2DArray().at(p).pixel( row, col )
			if badchs.pixel(row,col)>0:
				himg.SetBinContent(x+1,y+1,-50.0)
				continue			
			if val<10:
				continue
			himg.SetBinContent( x+1, y+1, val )
	himgs.append(himg)

# make path graphs
gpaths = []
for p in range(0,3):
	gpath = rt.TGraph( paths[p].size() )
	for i in range(0,paths[p].size()):
		node = paths[p].at(i)
		x = meta.pos_x( node.col )
		y = meta.pos_y( node.row )
		gpath.SetPoint( i, x, y )
	gpath.SetMarkerStyle(20)
	gpaths.append(gpath)


c = rt.TCanvas("c","c",1400,1200)
c.Divide(3,3)
c.Draw()

for p in range(0,3):
	c.cd( 3*p+1)
	himgs[p].Draw("COLZ")
	gpaths[p].Draw("LP")
	c.cd(3*p+2)
	hfscores[p].Draw("COLZ")
	c.cd(3*p+3)
	hgscores[p].Draw("COLZ")

c.Update()

raw_input()
