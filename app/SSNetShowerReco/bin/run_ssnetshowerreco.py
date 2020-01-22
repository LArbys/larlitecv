import os,sys,json,argparse

parser = argparse.ArgumentParser( description="Run shower reco and save to json file" )
parser.add_argument( "-ilcv", "--input-larcv",   type=str, required=True, help="Input larcv file. Should have ADC image, vertexer PGraph, SSNet images")
parser.add_argument( "-ill",  "--input-larlite", type=str, required=True, help="Input larlite file. Should have tracker trees")
parser.add_argument( "-o",    "--output-json",   type=str, required=True, help="Output JSON file" )

args = parser.parse_args()

import ROOT as rt
from ROOT import std
from larlite import larlite
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

data = {"entries":[]}

for ientry in xrange(nentries):
    print "[ENTRY ",ientry,"]"
    
    iolcv.read_entry(ientry)
    ioll.go_to(ientry)

    ok = showerreco.process( iolcv, ioll )

    entrydata = { "run":iolcv.event_id().run(),
                  "subrun":iolcv.event_id().subrun(),
                  "event":iolcv.event_id().event(),
                  "vertices":[] }
    
    for ivtx in xrange(showerreco.numVertices()):
        entrydata["vertices"].append( [ showerreco.getVertexShowerEnergy(ivtx,p) for p in xrange(3) ] )

    data["entries"].append( entrydata )

    if ientry>=2:
        print "break"
        break
    
print "output json"
fout = open(args.output_json, 'w' )
json.dump( data, fout )
fout.close()

print "close out"
iolcv.finalize()
ioll.close()


