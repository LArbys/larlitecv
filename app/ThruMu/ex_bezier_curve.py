import os,sys
import ROOT as rt
from ROOT import std
from larlitecv import larlitecv

c = rt.TCanvas("c1","c1",800,600)

cpts = std.vector("std::vector<float>")(4)

# PT1
cpts.at(0).resize(2,0)
cpts[0][0] = 1.0
cpts[0][1] = 1.0

# PT2
cpts.at(1).resize(2,0)
cpts[1][0] = 2.0
cpts[1][1] = 3.0

# PT3
cpts.at(2).resize(2,0)
cpts[2][0] = 4.0
cpts[2][1] = 3.0

# PT3
cpts.at(3).resize(2,0)
cpts[3][0] = 3.0
cpts[3][1] = 1.0

b = larlitecv.BezierCurve( cpts )

v = b.getCurve(100)

g = rt.TGraph(100)
for i in range(0,100):
    g.SetPoint(i,v.at(i).at(0), v.at(i).at(1))
box = rt.TGraph( 5 )
for i in range(0,4):
	box.SetPoint( i, cpts[i][0], cpts[i][1] )
box.SetPoint(4, cpts[0][0], cpts[0][1] )

c.Draw()
g.Draw("ALP")
box.SetLineColor(rt.kRed)
box.Draw("LP")
c.Update()
raw_input()



