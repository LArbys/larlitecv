import os,sys
import ROOT
from larcv import larcv
from larlitecv import larlitecv
from ROOT import std

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
valid_range.at(0).push_back(-100)
valid_range.at(0).push_back(1200)
valid_range.at(1).push_back(-150)
valid_range.at(1).push_back(150)

intersection3plane = std.vector("vector<int>")()
vertex3plane = std.vector("vector<float>")()
area3plane = std.vector("float")()
intersections2plane = std.vector("vector<int>")()
vertex2plane = std.vector("vector<float>")()


ubtool = larcv.UBWireTool.findWireIntersections( wirelists, valid_range, intersection3plane, vertex3plane, area3plane, intersections2plane, vertex2plane )

print "intersections: ",intersection3plane.size()
for i in range(0,intersection3plane.size()):
    print "intersection ",i,": ",intersection3plane.at(i).at(0), intersection3plane.at(i).at(1), intersection3plane.at(i).at(2)," area=",area3plane.at(i),
    print " vertex(z,y)=",vertex3plane.at(i).at(0),",",vertex3plane.at(i).at(1),","
    
