import ROOT as rt
from larlitecv import larlitecv 
from ROOT import std
import numpy as np
from array import array

rand = rt.TRandom3( 0 )


f = rt.TFile("temp.root", "RECREATE")
c = rt.TCanvas("c","c",1200,600)

h = rt.TH2D("h","h",100,0,1000,100,-120,120)

algo = larlitecv.FlashMuonTaggerAlgo(0)

while True:

    ls1 = std.vector("vector<float>")()
    start1 = std.vector("float")(2,0.0)
    start1[0] = np.random.uniform()*1000
    start1[1] = -115.0
    end1 = std.vector("float")(2,0.0)
    end1[0] = np.random.uniform()*1000
    end1[1] = +117.0
    ls1.push_back( start1 )
    ls1.push_back( end1 )

    ls2 = std.vector("vector<float>")()
    start2 = std.vector("float")(2,0.0)
    start2[0] = np.random.uniform()*1000
    start2[1] = -115.0
    end2 = std.vector("float")(2,0.0)
    end2[0] = np.random.uniform()*1000
    end2[1] = +117.0
    ls2.push_back( start2 )
    ls2.push_back( end2 )

    g1 = rt.TGraph(2)
    g1.SetPoint(0, start1[0], start1[1] )
    g1.SetPoint(1, end1[0], end1[1] )
    g2 = rt.TGraph(2)
    g2.SetPoint(0, start2[0], start2[1] )
    g2.SetPoint(1, end2[0], end2[1] )
    g2.SetLineColor(rt.kBlue)

    intersect = std.vector("float")(2,0.0)
    crosses = rt.Long(0)
    print crosses

    algo.lineSegmentIntersection2D( ls1, ls2, intersect, crosses )
    marker = rt.TMarker( intersect[0], intersect[1], 20 )
    marker.SetMarkerColor(rt.kRed)

    print crosses
    print " intersection point: ",intersect[0],intersect[1]

    c.Clear()
    h.Draw()
    g1.Draw("LP")
    g2.Draw("LP")
    marker.Draw("same")

    c.Update()


    print "[enter to continue]"
    raw_input()

