import os,sys
import ROOT
from ROOT import std
from larlite import larlite
from larlitecv import larlitecv

print "TEST RECLUSTER"

datacoord = larlitecv.DataCoordinator()
datacoord.add_inputfile( "../TaggerCROI/bin/tagger_anaout_larcv_seg.root",   "larcv" )
datacoord.add_inputfile( "../TaggerCROI/bin/tagger_anaout_larlite_seg.root", "larlite" )
datacoord.initialize()

nentries = datacoord.get_nentries("larlite")
print "NENTRIES: ", nentries

recluster_algo = larlitecv.Track3DRecluster()

for ientry in range(nentries):
    datacoord.goto_entry(ientry,"larlite")
    
    print "Entry: ",ientry
    ev_stopmu_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "stopmu3d" )
    ev_thrumu_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "thrumu3d" )

    print "StopMu Tracks: ",ev_stopmu_tracks.size()
    print "ThruMu Tracks: ",ev_thrumu_tracks.size()

    for itrack in range(ev_stopmu_tracks.size()+ev_thrumu_tracks.size()):
        if itrack<ev_stopmu_tracks.size():
            track = ev_stopmu_tracks[itrack]
        else:
            track = ev_thrumu_tracks[itrack-ev_stopmu_tracks.size()]
        npts = track.NumberTrajectoryPoints()
        #print "track %d number of pts="%(npts),npts
        path = std.vector("vector<float>")()
        for n in range(npts):
            pos_v = std.vector("float")(3)
            for i in range(3):
                pos_v[i] = track.LocationAtPoint(n)[i]
            #print "point ",n," pos_v=",pos_v[0],",",pos_v[1],",",pos_v[2]
            path.push_back( pos_v )
        print "track #%d path points: "%(itrack),path.size()
        recluster_algo.addPath( path )

    print "RECLUSTER"
    recluster_algo.recluster()
    
    break
