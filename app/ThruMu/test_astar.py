import os,sys,time
import numpy as np

print "Loading ROOT...",
s = time.time()
import ROOT as rt
from larlitecv import larlitecv
from larcv import larcv
print " ",time.time()-s,"secs"

rt.gStyle.SetOptStat(0)

def test_priority_queue():
    # test priority queue
    set = larlitecv.AStarSet()

    node10 = None
    for i in range(0,20):
        x = np.random.randint(0,100)
        y = np.random.randint(0,100)
        node = larlitecv.AStarNode( y, x )
        node.fscore = np.random.random()*100.0
        node.id = i
        print "push new node",node.id,node.col,node.row,node.fscore
        set.add( node )
        if i==10:
            node10 = node
        
    print "node10=",node10.id,node10.col,node10.row,node10.fscore
    print "CONTAINS node10: ",set.contains(node10)

            
    print "READ BACK"
    print "node10=",node10.id,node10.row,node10.col,node10.fscore
    #while not set.empty():
    #    node = set.pop()
    savenode = None
    for i in range(0,20):
        print "Check for node10 again: ",set.contains(node10)
        savenode = set.gettop()
        print savenode.id,savenode.row,savenode.col,savenode.fscore,savenode
        #set.pop() # get rid of the top
        if savenode==node10:
            print "node10"
            #print "popped node10: check again for it: ",set.contains(node10)
    print savenode


def test_astart_algo():
    # make a grid system

    grid = rt.TH2D( "grid","",500,0,500,500,0,500 )

    # define random circlular obstruction
    #r = 200.0
    r = 100.0
    x = np.random.random()*80 + 210.0
    y = np.random.random()*80 + 210.0

    print "fill all bins"

    meta = larcv.ImageMeta( 500, 500, 500, 500, 0, 500, 0 )
    img  = larcv.Image2D( meta )
    img.paint(0.0);
    s = time.time()
    for xbin in range(0,grid.GetXaxis().GetNbins()):
        for ybin in range(0,grid.GetYaxis().GetNbins()):
            if np.sqrt((xbin-x)*(xbin-x) + (ybin-y)*(ybin-y))>r:
                grid.SetBinContent( xbin+1, ybin+1, 10.0 )
                img.set_pixel( ybin, xbin, 10.0 )
    print "filled: ",time.time()-s,"secs"

    # algo
    config = larlitecv.AStarAlgoConfig()
    config.astar_threshold.push_back( 5.0 )
    config.astar_neighborhood.push_back( 2 )
    algo = larlitecv.AStarGridAlgo( config )
    algo.setVerbose(2)
    start = larlitecv.AStarNode(5,5)
    end   = larlitecv.AStarNode(495,495)

    grid.SetBinContent( start.col+1, start.row+1, 100.0 )
    grid.SetBinContent( end.col+1, end.row+1, 100.0 )

    print "find path ... "
    s = time.time()
    path = algo.findpath( img, 5, 5, 495, 495, 5.0 )
    print " done. ",time.time()-s,"secs"
    print "path points: ",path.size()
    for i in range(0,path.size()):
        node = path.at(i)
        grid.SetBinContent( node.col+1, node.row+1, 50.0 )

    c = rt.TCanvas("c","",800,600)
    c.Draw()
    grid.Draw("COLZ")
    c.Update()
    raw_input()

    

# test the AStarNode and AStarSet objects behave
#test_priority_queue()

# test astart algorithm
test_astart_algo()
