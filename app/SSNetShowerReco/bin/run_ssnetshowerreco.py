import os,sys,json,argparse
from math import sqrt

parser = argparse.ArgumentParser( description="Run shower reco and save to json file" )
parser.add_argument( "-ilcv", "--input-larcv",   type=str, required=True, help="Input larcv file. Should have ADC image, vertexer PGraph, SSNet images")
parser.add_argument( "-ill",  "--input-larlite", type=str, required=True, help="Input larlite file. Should have tracker trees")
parser.add_argument( "-o",    "--output-json",   type=str, required=True, help="Output JSON file" )
parser.add_argument( "-mc",   "--has-mc", default=False, action='store_true', help="Indicate input files have MC truth information" )

args = parser.parse_args()

import ROOT as rt
from ROOT import std
from larlite import larlite
from ROOT import larutil
larutil.SpaceChargeMicroBooNE
from larcv import larcv
from larlitecv import larlitecv


iolcv = larcv.IOManager( larcv.IOManager.kREAD, "lcvio" )
iolcv.add_in_file( args.input_larcv )
iolcv.initialize()

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename( args.input_larlite )
ioll.open()

nentries = iolcv.get_n_entries()

showerreco = larlitecv.ssnetshowerreco.SSNetShowerReco()
mcpg = larlitecv.mctruthtools.MCPixelPGraph()
sce  = larutil.SpaceChargeMicroBooNE() # larutil.SpaceChargeMicroBooNE.kMCC9_Forward )

data = {"entries":[]}

for ientry in xrange(nentries):
    print "[ENTRY ",ientry,"]"
    
    iolcv.read_entry(ientry)
    ioll.go_to(ientry)

    ok = showerreco.process( iolcv, ioll )

    entrydata = { "run":iolcv.event_id().run(),
                  "subrun":iolcv.event_id().subrun(),
                  "event":iolcv.event_id().event(),
                  "shower_energies":[],
                  "vertex_pos":[]}
    
    for ivtx in xrange(showerreco.numVertices()):
        entrydata["shower_energies"].append( [ showerreco.getVertexShowerEnergy(ivtx,p) for p in xrange(3) ] )
        entrydata["vertex_pos"].append( [ showerreco.getVertexPos(ivtx).at(p) for p in xrange(3) ] )

    # save vertex truth information
    if args.has_mc:
        # build graph and get primary particles
        mcpg.buildgraph( iolcv, ioll )
        mcpg.printGraph()
        node_v = mcpg.getPrimaryParticles()

        # determine topology, get electron energy if available
        nelectrons = 0
        nprotons   = 0
        nother     = 0
        pidX = []
        pidOther = []
        for inode in xrange(node_v.size()):
            node = node_v.at(inode)
            #print "node[",inode,"] pid=",node.pid," E=",node.E_MeV
            if abs(node.pid)==11:
                entrydata["true_electron_energy"] = node.E_MeV
                nelectrons += 1
            elif node.pid==2212:
                #print "found proton: ",node.E_MeV
                if node.E_MeV>60.0:
                    nprotons += 1
            else:
                if node.pid in [211,-211]:
                    # charged pions
                    nother += 1
                    pidOther.append( node.pid )
                else:
                    pidX.append(node.pid)
        entrydata["true_topology"] = "%de%dp%dX"%(nelectrons,nprotons,nother)
        print "topology",entrydata["true_topology"]
        #print "unknown PID: ",pidX

        # determine distance of reco vertices from true
        vtx_v = mcpg.findTrackID(-1).start
        entrydata["true_vertex"] = [ vtx_v[i] for i in xrange(3) ] # get the ROOT node
        offset_v = sce.GetPosOffsets( vtx_v[0], vtx_v[1], vtx_v[2] )
        vtx_sce_v = [ vtx_v[0]-offset_v[0]+0.7,
                      vtx_v[1]+offset_v[1],
                      vtx_v[2]+offset_v[2] ]
        entrydata["true_vertex_sce"] = vtx_sce_v

        entrydata["vertex_dist_from_truth"] = []
        for pos in entrydata["vertex_pos"]:
            d = 0.0
            for i in xrange(3):
                d += (pos[i]-vtx_sce_v[i])*(pos[i]-vtx_sce_v[i])
            d = sqrt(d)
            entrydata["vertex_dist_from_truth"].append(d)
                

    data["entries"].append( entrydata )

    #if ientry>=0:
    #    print "break"
    #    break
    
print "output json"
fout = open(args.output_json, 'w' )
json.dump( data, fout )
fout.close()

print "close out"
iolcv.finalize()
ioll.close()


