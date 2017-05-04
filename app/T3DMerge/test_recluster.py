import os,sys
import ROOT
from ROOT import std
from larlite import larlite
from larcv import larcv
from larlitecv import larlitecv

print "TEST RECLUSTER"

config = """Config: {
  IOManager: {
    IOMode: 2
    Verbosity: 2
    OutFileName: \"test_recluster_larcv.root\"
  }
  StorageManager: {
    IOMode: 2
    Verbosity: 2
    OutFileName: \"test_recluster_larlite.root\"
  }
}
"""

configout = """ConfigOut: {
  IOManager: {
    IOMode: 1
    Verbosity: 2
    OutFileName: \"out_recluster_larcv.root\"
  }
  StorageManager: {
    IOMode: 1
    Verbosity: 2
    OutFileName: \"out_recluster_larlite.root\"
  }
}
"""

fout = open("recluster.cfg",'w')
print >> fout,config
fout.close()

fout = open("out.cfg",'w')
print >> fout,configout
fout.close()

datacoord = larlitecv.DataCoordinator()
datacoord.configure( "recluster.cfg", "StorageManager", "IOManager", "Config" )
datacoord.add_inputfile( "../TaggerCROI/bin/tagger_anaout_larcv.root",   "larcv" )
datacoord.add_inputfile( "../TaggerCROI/bin/tagger_anaout_larlite.root", "larlite" )
datacoord.initialize()

dataout = larlitecv.DataCoordinator()
dataout.configure("out.cfg","StorageManager", "IOManager", "ConfigOut" )
dataout.initialize()

nentries = datacoord.get_nentries("larlite")
print "NENTRIES: ", nentries

recluster_algo = larlitecv.Track3DRecluster()
recluster_algo.setVerbosity(0)

thresholds = std.vector("float")(3,10.0)
neighborhood = std.vector("int")(3,5)
stepsize = 0.3

for ientry in range(0,5):
#for ientry in range(0,1):
    
    datacoord.goto_entry(ientry,"larlite")
    
    print "Entry: ",ientry
    ev_img_v         = datacoord.get_larcv_data( larcv.kProductImage2D, "modimg" )
    ev_gapch_v       = datacoord.get_larcv_data( larcv.kProductImage2D, "gapchs" )    
    ev_thrumu_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "thrumu3d" )    
    ev_stopmu_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "stopmu3d" )

    img_v   = ev_img_v.Image2DArray()
    gapch_v = ev_gapch_v.Image2DArray()
    
    print "StopMu Tracks: ",ev_stopmu_tracks.size()
    print "ThruMu Tracks: ",ev_thrumu_tracks.size()
    ntracks = ev_stopmu_tracks.size()+ev_thrumu_tracks.size()
    #ntracks = ev_thrumu_tracks.size()    
    print "NTRACKS=",ntracks
    
    for itrack in range(ntracks):
        if itrack<ev_thrumu_tracks.size():
            track = ev_thrumu_tracks[itrack]
        else:
            track = ev_stopmu_tracks[itrack-ev_thrumu_tracks.size()]
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
    recluster_tracks = recluster_algo.recluster()
    print "Number of recluster tracks: ",recluster_tracks.size()

    recluster_v = std.vector("larcv::Image2D")()
    for p in range(img_v.size()):
        recluster_img = larcv.Image2D( img_v[p].meta() )
        recluster_img.paint(0)
        recluster_v.push_back( recluster_img )
    print "Recluster images=",recluster_v.size()
    raw_input()    
    
    ev_recluster_out = dataout.get_larlite_data( larlite.data.kTrack, "recluster3d" )
    for itrack in range(recluster_tracks.size()):
        track = recluster_tracks[itrack]
        print "TRACK #%d"%(itrack),
        print " (",track.getPath().front()[0],",",track.getPath().front()[1],",",track.getPath().front()[2],") ",
        print " (",track.getPath().back()[0],",",track.getPath().back()[1],",",track.getPath().back()[2],") ",
        print " size=",track.getPath().size()
        #for ipt in range(track.getPath().size()):
        #    print "   [ipt] (",track.getPath()[ipt][0],",",track.getPath()[ipt][1],",",track.getPath()[ipt][2],")"
        lltrack = larlitecv.T3D2LarliteTrack( track )

        pix_clusters_v = track.getPixelsFromImages( img_v, gapch_v, thresholds, neighborhood, stepsize )
        for p in range(recluster_v.size()):
            for ipix in range(pix_clusters_v[p].size()):
                pix = pix_clusters_v[p][ipix]
                recluster_v[p].set_pixel( pix.Y(), pix.X(), 3*itrack )
        
        ev_recluster_out.push_back( lltrack )

    ev_recluster_pixels = dataout.get_larcv_data( larcv.kProductImage2D, "reclusterpixels" )    
    for p in range(recluster_v.size()):
        ev_recluster_pixels.Append( recluster_v[p] )
    print "Saving ",ev_recluster_pixels.Image2DArray().size()," recluster images"

    run    = datacoord.run()
    subrun = datacoord.subrun()
    event  = datacoord.event()
    dataout.set_id( run, subrun, event );
    dataout.save_entry()
    datacoord.save_entry()    
    break

dataout.finalize()
datacoord.finalize()
