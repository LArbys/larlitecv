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

start_tick = 2412
start_wire = 60
end_tick   = 3342
end_wire   = 10
use_bad_chs = False

config  = larlitecv.AStarDirAlgoConfig()
config.astar_threshold.resize(3,10.0)
config.astar_neighborhood.resize(3,5)
config.astar_start_padding = 2
config.astar_end_padding = 2
config.image_padding = 20

algo = larlitecv.AStarDirAlgo( config )
algo.setVerbose(0)

meta = ev_img.Image2DArray().at(2).meta()
badchs = ev_badchs.Image2DArray().at(2)
algo.setBadChImage( ev_badchs.Image2DArray().at(2) )
path = algo.findpath( ev_img.Image2DArray().at(2), meta.row(start_tick), meta.col(start_wire), meta.row(end_tick), meta.col(end_wire), 10.0, use_bad_chs )
fimg = algo.getScoreImages().at(0)
gimg = algo.getScoreImages().at(1)

print "path size: ",path.size()

# fill fscore image
print "fimg: ",fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y()
hfscore = rt.TH2D("fscores","fscores", fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y() )
hgscore = rt.TH2D("gscores","gscores", fimg.meta().cols(), fimg.meta().min_x(), fimg.meta().max_x(), fimg.meta().rows(), fimg.meta().min_y(), fimg.meta().max_y() )
for r in range(0,fimg.meta().rows()):
	for c in range(0,fimg.meta().cols()):
		rbin = fimg.meta().row( hfscore.GetYaxis().GetBinLowEdge(r+1)+0.1 )

		hfscore.SetBinContent( c+1, r+1, fimg.pixel(rbin,c) )
		hgscore.SetBinContent( c+1, r+1, gimg.pixel(rbin,c) )		

c = rt.TCanvas("c","c",1400,1200)
c.Divide(2,2)

low_wire = start_wire
high_wire = end_wire
if low_wire>high_wire:
	low_wire = end_wire
	high_wire = start_wire
high_wire += config.image_padding*meta.pixel_width()
low_wire  -= config.image_padding*meta.pixel_width()
if low_wire<0:
	low_wire = 0
print "wire range: ",low_wire,high_wire

nxbins = int( fabs(high_wire-low_wire)/meta.pixel_width() )
nybins = int( (end_tick-start_tick)/meta.pixel_height() )

print "bins x=",nxbins," y=",nybins
print "endy=",start_tick + meta.pixel_height()*nybins
print "endx=",low_wire + meta.pixel_width()*nxbins

himg = rt.TH2D("himg","himg", nxbins, low_wire, low_wire+meta.pixel_width()*nxbins, nybins, start_tick, start_tick + meta.pixel_height()*nybins )

for x in range(0,nxbins):
	for y in range(0,nybins):
		col = meta.col( himg.GetXaxis().GetBinLowEdge(x+1) )
		row = meta.row( himg.GetYaxis().GetBinLowEdge(y+1) )
		val = ev_img.Image2DArray().at(2).pixel( row, col )
		himg.SetBinContent( x+1, y+1, val )
		if badchs.pixel(row,col)>0:
			himg.SetBinContent(x+1,y+1,-10.0)

gpath = rt.TGraph( path.size() )
for i in range(0,path.size()):
	node = path.at(i)
	x = meta.pos_x( node.col )
	y = meta.pos_y( node.row )
	print i,x,y
	gpath.SetPoint( i, x, y )
gpath.SetMarkerStyle(20)

c.Draw()
c.cd(1)
himg.Draw("COLZ")
gpath.Draw("LP")
c.cd(2)
hfscore.Draw("COLZ")
c.cd(3)
hgscore.Draw("COLZ")

c.Update()

raw_input()