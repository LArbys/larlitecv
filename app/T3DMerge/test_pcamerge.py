import os,sys
import ROOT
from ROOT import std
from larlite import larlite
from larcv import larcv
from larlitecv import larlitecv
from math import sqrt

print "TEST RECLUSTER"

config = """Config: {
  IOManager: {
    IOMode: 2
    Verbosity: 2
    OutFileName: \"out_pcamerged_larcv.root\"
  }
  StorageManager: {
    IOMode: 2
    Verbosity: 2
    OutFileName: \"out_pcamerged_larlite.root\"
  }
}
"""

configout = """ConfigOut: {
  IOManager: {
    IOMode: 1
    Verbosity: 2
    OutFileName: \"out_pcamerged2_larcv.root\"
  }
  StorageManager: {
    IOMode: 1
    Verbosity: 2
    OutFileName: \"out_pcamerged2_larlite.root\"
  }
}
"""

fout = open("pcamerge.cfg",'w')
print >> fout,config
fout.close()

fout = open("pcamergeout.cfg",'w')
print >> fout,configout
fout.close()

datacoord = larlitecv.DataCoordinator()
datacoord.configure( "pcamerge.cfg", "StorageManager", "IOManager", "Config" )
datacoord.add_inputfile( "testfiles/tagger_anaout_larcv_mergetest.root",   "larcv" )
datacoord.add_inputfile( "testfiles/tagger_anaout_larlite_mergetest.root", "larlite" )
datacoord.initialize()

dataout = larlitecv.DataCoordinator()
dataout.configure("pcamergeout.cfg","StorageManager", "IOManager", "ConfigOut" )
dataout.initialize()

nentries = datacoord.get_nentries("larlite")
print "NENTRIES: ", nentries

recluster_algo = larlitecv.Track3DRecluster()
recluster_algo.setVerbosity(0)

recluster_algo2 = larlitecv.Track3DRecluster()
recluster_algo2.setVerbosity(0)

pcamerge_algo = larlitecv.T3DPCMerge()

thresholds = std.vector("float")(3,10.0)
neighborhood = std.vector("int")(3,5)
stepsize = 0.3

for ientry in range(0,1):
    
    datacoord.goto_entry(ientry,"larlite")
    
    print "Entry: ",ientry
    ev_img_v           = datacoord.get_larcv_data( larcv.kProductImage2D, "modimg" )
    ev_gapch_v         = datacoord.get_larcv_data( larcv.kProductImage2D, "gapchs" )    
    ev_thrumu_tracks   = datacoord.get_larlite_data( larlite.data.kTrack, "thrumu3d" )    
    ev_stopmu_tracks   = datacoord.get_larlite_data( larlite.data.kTrack, "stopmu3d" )
    ev_untagged_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "untagged3d" )
    ev_streclustered_tracks = datacoord.get_larlite_data( larlite.data.kTrack, "streclustered3d" )    

    img_v   = ev_img_v.Image2DArray()
    gapch_v = ev_gapch_v.Image2DArray()
    
    print "StopMu Tracks: ",ev_stopmu_tracks.size()
    print "ThruMu Tracks: ",ev_thrumu_tracks.size()
    print "Untagged Tracks: ",ev_untagged_tracks.size()
    print "Stop/Thru Reclustered Tracks",ev_streclustered_tracks.size()

    # recluster untagged
    for itrack in range(ev_untagged_tracks.size()):
        track = ev_untagged_tracks[itrack]
        npts = track.NumberTrajectoryPoints()
        path = std.vector("vector<double>")()
        for n in range(npts):
            pos_v = std.vector("double")(3)
            for i in range(3):
                pos_v[i] = track.LocationAtPoint(n)[i]
            path.push_back( pos_v )
        recluster_algo.addPath( path )


    reclustered_untagged = recluster_algo.recluster()
    print "RECLUSTER UNTAGGED: ",reclustered_untagged.size()

    for trackset in [ reclustered_untagged, ev_streclustered_tracks ]:
        ntracks = trackset.size()
        if trackset==reclustered_untagged:
            for itrack in range(ntracks):
                track = trackset[itrack]                    
                #recluster_algo.addTrack(  track )
                print "track #%d path points: "%(itrack),track.getPath().size()                
                recluster_algo2.addPath( track.getPath() )
        else:
            for itrack in range(ntracks):
                track = trackset[itrack]                                    
                npts = track.NumberTrajectoryPoints()
                #print "track %d number of pts="%(npts),npts
                path = std.vector("vector<double>")()
                for n in range(npts):
                    pos_v = std.vector("double")(3)
                    for i in range(3):
                        pos_v[i] = track.LocationAtPoint(n)[i]
                    #print "point ",n," pos_v=",pos_v[0],",",pos_v[1],",",pos_v[2]
                    path.push_back( pos_v )
                print "track #%d path points: "%(itrack),path.size()
                #t3d = larlitecv.T3DCluster( path )
                #recluster_algo.addTrack( t3d )
                recluster_algo2.addPath( path )

    t3d_v = recluster_algo2.recluster()
    print "RECLUSTER UNTAGGED + STRECLUSTERED: ",t3d_v.size()
    
    # t3d_v = std.vector("larlitecv::T3DCluster")()
    # for trackset in [ reclustered_untagged, ev_streclustered_tracks ]:
    #     ntracks = trackset.size()
    #     if trackset==reclustered_untagged:
    #         for itrack in range(ntracks):
    #             track = trackset[itrack]                    
    #             t3d_v.push_back( track )
    #     else:
    #         for itrack in range(ntracks):
    #             track = trackset[itrack]                                    
    #             npts = track.NumberTrajectoryPoints()
    #             #print "track %d number of pts="%(npts),npts
    #             path = std.vector("vector<double>")()
    #             for n in range(npts):
    #                 pos_v = std.vector("double")(3)
    #                 for i in range(3):
    #                     pos_v[i] = track.LocationAtPoint(n)[i]
    #                 #print "point ",n," pos_v=",pos_v[0],",",pos_v[1],",",pos_v[2]
    #                 path.push_back( pos_v )
    #             print "track #%d path points: "%(itrack),path.size()
    #             t3d = larlitecv.T3DCluster( path )
    #             t3d_v.push_back( t3d )

    print "TOTAL LOADED TRACKS: ",t3d_v.size()
    
    # Check PCA Track
    #pcadist = 0
    #for i in range(3):
    #    pcadist += t3d_v.front().getPCADir(0)[i]*t3d_v.front().getPCADir(0)[i]
    #pcadist = sqrt(pcadist)
    #print "test norm: ",pcadist
    #print "test track: (",t3d_v.front().getPCADir(0)[0],t3d_v.front().getPCADir(0)[1],t3d_v.front().getPCADir(0)[2],")"
    #print "start: (",t3d_v.front().getPath().front()[0],t3d_v.front().getPath().front()[1],t3d_v.front().getPath().front()[2],")"
    #print "end: (",t3d_v.front().getPath().back()[0],t3d_v.front().getPath().back()[1],t3d_v.front().getPath().back()[2],")"

    pcmerged_v = pcamerge_algo.merge( t3d_v )
    #pcmerged_v = t3d_v
    print "check merged: ",pcmerged_v.size(), " vs ", t3d_v.size()
    
    merged_v = std.vector("larcv::Image2D")()
    for p in range(img_v.size()):
        merged_img = larcv.Image2D( img_v[p].meta() )
        merged_img.paint(0)
        merged_v.push_back( merged_img )
    print "Merged images=",merged_v.size()
    raw_input()    
    
    ev_merged_out = dataout.get_larlite_data( larlite.data.kTrack, "pcamerged3d" )
    for itrack in range(pcmerged_v.size()):
        track = pcmerged_v[itrack]
        print "TRACK #%d"%(itrack),
        print " (",track.getPath().front()[0],",",track.getPath().front()[1],",",track.getPath().front()[2],") ",
        print " (",track.getPath().back()[0],",",track.getPath().back()[1],",",track.getPath().back()[2],") ",
        print " size=",track.getPath().size()
        #for ipt in range(track.getPath().size()):
        #    print "   [ipt] (",track.getPath()[ipt][0],",",track.getPath()[ipt][1],",",track.getPath()[ipt][2],")"
        lltrack = larlitecv.T3D2LarliteTrack( track )

        pix_clusters_v = track.getPixelsFromImages( img_v, gapch_v, thresholds, neighborhood, stepsize )
        for p in range(merged_v.size()):
            for ipix in range(pix_clusters_v[p].size()):
                pix = pix_clusters_v[p][ipix]
                merged_v[p].set_pixel( pix.Y(), pix.X(), 3*itrack )
        
        ev_merged_out.push_back( lltrack )

    ev_merged_pixels = dataout.get_larcv_data( larcv.kProductImage2D, "pcamergedpixels" )    
    for p in range(merged_v.size()):
        ev_merged_pixels.Append( merged_v[p] )
    print "Saving ",ev_merged_pixels.Image2DArray().size()," merged images"

    run    = datacoord.run()
    subrun = datacoord.subrun()
    event  = datacoord.event()
    dataout.set_id( run, subrun, event );
    dataout.save_entry()
    datacoord.save_entry()    
    break

dataout.finalize()
datacoord.finalize()
