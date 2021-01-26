#!/usr/bin/env python
import os,sys,json,argparse
from math import sqrt

parser = argparse.ArgumentParser( description="Run shower reco and save to json file" )
parser.add_argument( "-ill_list",   "--input-larlite", type=str,  default=None,   help="Input larlite file list. Should have MC.")
parser.add_argument( "-o",     "--output",        type=str,  required=True,  help="Output file name. if both, the stem for both." )
args = parser.parse_args()

import ROOT as rt
from ROOT import std
from larlite import larlite
from ROOT import larutil
larutil.SpaceChargeMicroBooNE
from larcv import larcv
from larlitecv import larlitecv

data = {"entries":[]}
# json output
jout_name = args.output
jout_name += ".json"



ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
if args.input_larlite is not None:
    ioll.add_in_filename(args.input_larlite)
    ioll.open()

# larcv for r,s,e info
iolcv = larcv.IOManager( larcv.IOManager.kREAD, "lcvio" )
iolcv.add_in_file( args.input_larlite)
iolcv.initialize()


GetTruthVariables = larlitecv.ssnetshowerreco.GetTruthVariables()
GetTruthVariables.initialize()


nentries = ioll.get_entries()
for ientry in xrange(nentries):
    print "[ENTRY ",ientry,"]"

    ioll.go_to(ientry)
    iolcv.read_entry(ientry)

    ok = GetTruthVariables.process( iolcv,ioll, ientry )
    print("larlite run", ioll.run_id(), "larcv run", iolcv.event_id().run())
    print("larlite subrun", ioll.subrun_id(), "larcv subrun", iolcv.event_id().subrun())
    print("larlite event", ioll.event_id(), "larcv eventrun", iolcv.event_id().event())
    ev_wire = iolcv.get_data(larcv.kProductImage2D,"wire")
    print(ev_wire.run(),ev_wire.subrun(), ev_wire.event())

    # print (iolcv.event_id().event())
    entrydata = { "run":ioll.run_id(),
                  "subrun":ioll.subrun_id(),
                  "event":ioll.event_id(),
                  "SSNetShowerAverage":GetTruthVariables.getSSNetShowerAverage(),
                  "particle_pdg":[],
                  "particle_status":[],
                  "particle_Energy":[],
                  "particle_Px":[],
                  "particle_Py":[],
                  "particle_Pz":[],
                  }

    # Save first shower output
    for i in xrange(GetTruthVariables.NumParts()):
        entrydata["particle_pdg"].append( GetTruthVariables.getParticlePDG(i) )
        entrydata["particle_status"].append( GetTruthVariables.getParticlePDG(i) )
        entrydata["particle_Energy"].append( GetTruthVariables.getPartE(i) )
        entrydata["particle_Px"].append( GetTruthVariables.getPartMomentumX(i) )
        entrydata["particle_Py"].append( GetTruthVariables.getPartMomentumY(i) )
        entrydata["particle_Pz"].append( GetTruthVariables.getPartMomentumZ(i) )


    data["entries"].append( entrydata )

print "output json"

GetTruthVariables.finalize()

print "close out"

ioll.close()
iolcv.finalize()

fout = open(jout_name, 'w' )
json.dump( data, fout )
fout.close()
