import os,sys
import ROOT
from larcv import larcv
from larlitecv import larlitecv
from ROOT import std
from ROOT import Long, Double

if len(sys.argv)==4:
    u = int(sys.argv[1])
    v = int(sys.argv[2])
    y = int(sys.argv[3])
else:
    u = 805
    v = 1189
    y = 1322

wirelists = std.vector("vector<int>")(3)
wirelists.at(0).push_back(u)
wirelists.at(1).push_back(v)
wirelists.at(2).push_back(y)

valid_range = std.vector("vector<float>")(2)
valid_range.at(0).push_back(-500)
valid_range.at(0).push_back(1500)
valid_range.at(1).push_back(-200)
valid_range.at(1).push_back(200)

intersection3plane = std.vector("vector<int>")()
vertex3plane = std.vector("vector<float>")()
area3plane = std.vector("float")()
intersections2plane = std.vector("vector<int>")()
vertex2plane = std.vector("vector<float>")()

# old way
#ubtool = larcv.UBWireTool.findWireIntersections( wirelists, valid_range, intersection3plane, vertex3plane, area3plane, intersections2plane, vertex2plane )

wids = std.vector("int")(3,0)
wids[0] = u
wids[1] = v
wids[2] = y

poszy = std.vector("float")(2,0.0)
crosses = Long(0)
print wids,poszy
area3plane.push_back(0.0)
area = Double(0)

larcv.UBWireTool.wireIntersection( wids, poszy, area, crosses )

print "Crosses: ",crosses
print "intersection : (z,y)=",poszy.at(0), poszy.at(1)," area=",area
    
