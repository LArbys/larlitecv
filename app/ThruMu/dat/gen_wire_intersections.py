import os,sys,time
s = time.time()
try:
   import cPickle as pickle
except:
   import pickle
from math import fabs
import ROOT
from larlite import larlite
from larcv import larcv
from ROOT import larutil
from ROOT import std
print "loading modules: ",time.time()-s,"secs"

print "GEN WIRE INTERSECTIONS"
ydata = larcv.UBWireTool.getWireData(2)
udata = larcv.UBWireTool.getWireData(0)
vdata = larcv.UBWireTool.getWireData(1)

mode = "tight"

# Top matches
# Bottom matches
top_matches = []
bot_matches = []
top_positions = []
bot_positions = []
for yid in range(0,3456):
    #print "Y-wire=",yid
    pos = std.vector("double")(3,0.0)
    pos[2] = ydata.wireStart[yid].at(2)

    # Top
    if mode=="tight":
       pos[1] = 104.0 # tight
    else:
       pos[1]  = 114.0 # loose
    top_uwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 0 ))
    top_vwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 1 ))
    top_match = ( int(top_uwire), int(top_vwire), yid)
    top_matches.append( top_match )
    top_positions.append( (pos[0],pos[1],pos[2]) )
    
    # bottom
    if mode=="tight":
       pos[1] = -100.0 # tight
    else:
       pos[1] = -114.0 # loose
    bot_uwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 0 ))
    bot_vwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 1 ))
    bot_match = ( int(bot_uwire), int(bot_vwire), yid)
    bot_matches.append( bot_match )
    bot_positions.append( (pos[0],pos[1],pos[2]) )

    #print "  top: ",top_match
    #print "  bottom: ",bot_match

# Upstream
upstream_matches = []
upstream_positions = []
for uid in range(0,2400):
    ustart = udata.wireStart[uid]
    if ustart[2]>0.036+10.0:
        continue
    udir   = udata.wireDir[uid]
    if mode=="tight":
       z = 0.036+10.0
    else:
       z = 0.036+5.0
    y = ustart[1]+udir[1]*fabs(z-ustart[2])
    pos = std.vector("double")(3,0.0)
    pos[1] = y
    pos[2] = z
    up_vwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 1 ))
    up_ywire = round(larutil.Geometry.GetME().WireCoordinate( pos, 2 ))
    
    up_match = (uid,int(up_vwire),int(up_ywire))
    print "Upstream match: uid=",uid," ",up_match," z=",ustart[2]
    upstream_matches.append( up_match )
    upstream_positions.append( (pos[0],pos[1],pos[2]) )

# Downstream
downstream_matches = []
downstream_positions = []
for vid in range(0,2400):
    vstart = vdata.wireStart[vid]
    if vstart[2]<1036.96472168-10.0:
        continue
    vdir   = vdata.wireDir[vid]
    if mode=="tight":
       z = 1036.96472168-10.0
    else:
       z = 1036.96472168-5.0
    y = vstart[1]+vdir[1]*fabs(vstart[2]-z)
    pos = std.vector("double")(3,0.0)
    pos[1] = y
    pos[2] = z
    down_uwire = round(larutil.Geometry.GetME().WireCoordinate( pos, 0 ))
    down_ywire = round(larutil.Geometry.GetME().WireCoordinate( pos, 2 ))
    
    down_match = (int(down_uwire),vid,int(down_ywire))
    print "Downstream match: vid=",vid," ",down_match," z=",vstart[2]
    downstream_matches.append( down_match )
    downstream_positions.append( (pos[0],pos[1],pos[2]) )

for matches in [("top",top_matches),("bottom",bot_matches),("upstream",upstream_matches),("downstream",downstream_matches)]:
   f = open( "%s_matches_sce_%s.pickle"%(matches[0],mode), 'w' )
   pickle.dump( matches[1], f )
   f.close()

for positions in [("top",top_positions),("bottom",bot_positions),("upstream",upstream_positions),("downstream",downstream_positions)]:
   f = open( "%s_positions_sce_%s.pickle"%(positions[0],mode), 'w' )
   pickle.dump( positions[1], f )
   f.close()
