import os,sys,json,argparse
from math import sqrt

parser = argparse.ArgumentParser( description="Run shower reco and save to json file" )
parser.add_argument( "-ilcv",  "--input-larcv",   type=str,  required=True,  help="Input larcv file. Should have ADC image, vertexer PGraph, SSNet images")
parser.add_argument( "-ilcvtruth",  "--input-larcvtruth",   type=str,  required=False,  help="Input larcv file. Should have segment images")
parser.add_argument( "-iimgs", "--input-images",  type=str,  default=None,   help="Input image file. Should be a dlmerged file for uncalibrated images or a calibrated* file in stage1")
parser.add_argument( "-ill",   "--input-larlite", type=str,  default=None,   help="Input larlite file. Should have tracker trees and MC.")
parser.add_argument( "-f",     "--output-format", type=str,  default="json", help="Set output format. Options={'json','larlite','both'}")
parser.add_argument( "-o",     "--output",        type=str,  required=True,  help="Output file name. if both, the stem for both." )
parser.add_argument( "-adc",   "--adc-tree",      type=str,  default="wire", help="Name of tree containing ADC images. [ Default: 'wire' ]")
parser.add_argument( "-cal",   "--use-calib", default=False, action='store_true', help="Use calibrated conversions")
parser.add_argument( "-mc",    "--has-mc", default=False, action='store_true', help="Indicate input files have MC truth information" )
parser.add_argument( "-sec",   "--second-shower", default=False, action='store_true', help="Run second shower search")
parser.add_argument( "-ncpi0",   "--use-ncpi0", default=False, action='store_true', help="Using NCPi0 true info")
parser.add_argument( "-nueint",   "--use-nueint", default=False, action='store_true', help="Using Nueint true info")
args = parser.parse_args()

output_formats = ['json','larlite','both']
if args.output_format not in output_formats:
    raise ValueError("invalid output format specified. options: {}".format(output_formats))

import ROOT as rt
from ROOT import std
from larlite import larlite
from ROOT import larutil
larutil.SpaceChargeMicroBooNE
from larcv import larcv
from larlitecv import larlitecv


iolcv = larcv.IOManager( larcv.IOManager.kREAD, "lcvio" )
iolcv.add_in_file( args.input_larcv )
if args.input_larcvtruth is not None:
    iolcv.add_in_file( args.input_larcvtruth)
if args.input_images is not None:
    ilcv.add_in_file( args.input_images )
iolcv.initialize()

# deprecated
#ioimgs = larcv.IOManager( larcv.IOManager.kREAD, "imgsio" )
#ioimgs.add_in_file( args.input_images )
#ioimgs.initialize()

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
if args.input_larlite is not None:
    ioll.add_in_filename( args.input_larlite )
    ioll.open()

# larlite output
llout_name = args.output
if args.output_format in ['larlite','both']:

    if args.output_format=='both':
        llout_name += "_larlite.root"

    outll = larlite.storage_manager( larlite.storage_manager.kWRITE )
    outll.set_out_filename( llout_name )
    outll.open()
else:
    outll = None

# json output
jout_name = args.output
if args.output_format=='both':
    jout_name += ".json"

nentries = iolcv.get_n_entries()

showerreco = larlitecv.ssnetshowerreco.SSNetShowerReco()
mcpg = larlitecv.mctruthtools.MCPixelPGraph()
sce  = larutil.SpaceChargeMicroBooNE() # larutil.SpaceChargeMicroBooNE.kMCC9_Forward

showerreco.set_adc_treename( args.adc_tree )
if args.use_calib:
    showerreco.use_calibrated_pixsum2mev( True )
if args.second_shower:
    showerreco.use_second_shower( True )
if args.use_ncpi0:
    showerreco.use_ncpi0( True )
if args.use_nueint:
    showerreco.use_nueint( True )

data = {"entries":[]}

for ientry in xrange(nentries):
    print "[ENTRY ",ientry,"]"

    iolcv.read_entry(ientry)
    if args.input_larlite is not None:
        ioll.go_to(ientry)

    ok = showerreco.process( iolcv, ioll, ientry )
    if outll is not None:
        showerreco.store_in_larlite( outll )
        outll.set_id( iolcv.event_id().run(), iolcv.event_id().subrun(), iolcv.event_id().event() )
        outll.next_event()

    print (iolcv.event_id().event())
    entrydata = { "run":iolcv.event_id().run(),
                  "subrun":iolcv.event_id().subrun(),
                  "event":iolcv.event_id().event(),
                  "shower_energies":[],
                  "shower_sumQs":[],
                  "shower_shlengths":[],
                  "vertex_pos":[]}

    for ivtx in xrange(showerreco.numVertices()):
        entrydata["shower_energies"].append( [ showerreco.getVertexShowerEnergy(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_sumQs"].append( [ showerreco.getVertexShowerSumQ(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_shlengths"].append( [ showerreco.getVertexShowerShlength(ivtx,p) for p in xrange(3) ] )
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

    if args.second_shower:
        entrydata["secondshower_energies"] = []
        entrydata["secondshower_sumQs"] = []
        entrydata["secondshower_shlengths"] = []

        for ivtx in xrange(showerreco.numVertices()):
            entrydata["secondshower_energies"].append( [ showerreco.getVertexSecondShowerEnergy(ivtx,p) for p in xrange(3) ] )
            entrydata["secondshower_sumQs"].append( [ showerreco.getVertexSecondShowerSumQ(ivtx,p) for p in xrange(3) ] )
            entrydata["secondshower_shlengths"].append( [ showerreco.getVertexSecondShowerShlength(ivtx,p) for p in xrange(3) ] )

    if args.use_ncpi0:
        entrydata["true_shower_energies"] = []
        entrydata["true_shower_starts"] = []
        entrydata["remaining_adc"] = []
        entrydata["overlap_fraction1"] = []
        entrydata["overlap_fraction2"] = []
        entrydata["truth_match"] = []




        for ivtx in xrange(showerreco.numShowers()):
            entrydata["truth_match"].append( [showerreco.getShowerTruthMatch(shower) for shower in xrange(6)])
            entrydata["true_shower_energies"].append( [ showerreco.getTrueShowerEnergy(ivtx) for shower in xrange(2) ] )
            entrydata["true_shower_starts"].append( [ showerreco.getTrueShowerStarts(ivtx).at(p) for p in xrange(3) ] )
            entrydata["remaining_adc"].append( [showerreco.getRemainingADC()])
            entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(plane,0) for plane in xrange(2) ] )
            entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(plane,1) for plane in xrange(2) ] )
            entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(plane,0) for plane in xrange(2) ] )
            entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(plane,1) for plane in xrange(2) ] )

    if args.use_nueint:
        entrydata["uplane_profile"] = []
        entrydata["vplane_profile"] = []
        entrydata["yplane_profile"] = []

        for ii in xrange(showerreco.numpointsU()):
            entrydata["uplane_profile"].append( [ showerreco.getUPlaneShowerProfile(ii,index) for index in xrange(2) ] )
        for ii in xrange(showerreco.numpointsV()):
            entrydata["vplane_profile"].append( [ showerreco.getVPlaneShowerProfile(ii,index) for index in xrange(2) ] )
        for ii in xrange(showerreco.numpointsY()):
            entrydata["yplane_profile"].append( [ showerreco.getYPlaneShowerProfile(ii,index) for index in xrange(2) ] )

    #if ientry>=0:
    #    print "break"
    #    break

    data["entries"].append( entrydata )

print "output json"
if args.output_format in ['json','both']:
    fout = open(jout_name, 'w' )
    json.dump( data, fout )
    fout.close()

print "close out"
outll.close()
iolcv.finalize()
if args.input_larlite is not None:
    ioll.close()
