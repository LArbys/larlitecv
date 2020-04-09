#!/usr/bin/env python
import os,sys,json,argparse
from math import sqrt

parser = argparse.ArgumentParser( description="Run shower reco and save to json file" )
parser.add_argument( "-ilcv",  "--input-larcv",   type=str,  required=True,  help="Input larcv file. Should have ADC image, vertexer PGraph, SSNet images [REQUIRED]")
parser.add_argument( "-ilcvtruth",  "--input-larcvtruth",   type=str,  required=False,  help="Input larcv file. Should have segment images [OPTIONAL]")
parser.add_argument( "-iimgs", "--input-images",  type=str,  default=None,   help="Input image file. Should be a dlmerged file for uncalibrated images or a calibrated* file in stage1")
parser.add_argument( "-ill",   "--input-larlite", type=str,  default=None,   help="Input larlite file. Should have tracker trees and MC.")
parser.add_argument( "-f",     "--output-format", type=str,  default="json", help="Set output format. Options={'json','larlite','both'}")
parser.add_argument( "-o",     "--output",        type=str,  required=True,  help="Output file name. if both, the stem for both." )
parser.add_argument( "-adc",   "--adc-tree",      type=str,  default="wire", help="Name of tree containing ADC images. [ Default: 'wire' ]")
parser.add_argument( "-ssnet",   "--ssnet-tree",      type=str,  default="uburn", help="Name of tree containing SSNet images. [ Default: 'uburn' ]")
parser.add_argument( "-cal",   "--use-calib", default=False, action='store_true', help="Use calibrated conversions")
parser.add_argument( "-mc",    "--has-mc", default=False, action='store_true', help="Indicate input files have MC truth information" )
parser.add_argument( "-sec",   "--second-shower", default=False, action='store_true', help="Run second shower search")
parser.add_argument( "-ncpi0",   "--use-ncpi0", default=False, action='store_true', help="Using NCPi0 true info")
parser.add_argument( "-nueint",   "--use-nueint", default=False, action='store_true', help="Using Nueint true info")
parser.add_argument( "-bnb",   "--use-bnb", default=False, action='store_true', help="Using Bnb Overlay info")
parser.add_argument( "-tn",    "--tree-name", default="ssnetshowerreco", type=str, help="Name of the output trees [ Default: 'ssnetshowerreco' ]")
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
    iolcv.add_in_file( args.input_images )
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

uselarlite =False
if args.output_format in ['larlite','both']:
    uselarlite =True

showerreco = larlitecv.ssnetshowerreco.SSNetShowerReco(uselarlite,llout_name)
showerreco.set_adc_treename( args.adc_tree )
showerreco.set_ssnet_shower_stemname( args.ssnet_tree )
if args.use_calib:
    print('USING CALIB')
    showerreco.use_calibrated_pixsum2mev( True )
if args.second_shower:
    showerreco.use_second_shower( True )
if args.use_ncpi0:
    showerreco.use_ncpi0( True )
if args.use_nueint:
    showerreco.use_nueint( True )
if args.use_bnb:
    showerreco.use_bnb( True )
if args.has_mc:
    showerreco.use_mc( True )

showerreco.set_output_treename( args.tree_name )
showerreco.initialize()

mcpg = larlitecv.mctruthtools.MCPixelPGraph()
sce  = larutil.SpaceChargeMicroBooNE() # larutil.SpaceChargeMicroBooNE.kMCC9_Forward

data = {"entries":[]}

nentries = iolcv.get_n_entries()
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

    # print (iolcv.event_id().event())
    entrydata = { "run":iolcv.event_id().run(),
                  "subrun":iolcv.event_id().subrun(),
                  "event":iolcv.event_id().event(),
                  "shower_energies":[],
                  "shower_sumQs":[],
                  "shower_fullSumQs":[],
                  "shower_smallSumQs":[],
                  "shower_triAreas":[],
                  "shower_shlengths":[],
                  "vertex_pos":[],
                  "shower_gap":[],
                  "shower_direction_3d":[],
                  "shower_direction_2d":[],
                  "shower_opening_2d":[],
                  "shower_start_2d":[],
                  }

    # Save first shower output
    for ivtx in xrange(showerreco.numVertices()):
        entrydata["shower_energies"].append( [ showerreco.getVertexShowerEnergy(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_sumQs"].append( [ showerreco.getVertexShowerSumQ(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_fullSumQs"].append( [ showerreco.getVertexShowerFullSumQ(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_smallSumQs"].append( [ showerreco.getVertexShowerSmallSumQ(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_triAreas"].append( [ showerreco.getVertexShowerTriArea(ivtx,p) for p in xrange(3) ] )
        for p in xrange(3):
          print 'Python full sum:',showerreco.getVertexShowerFullSumQ(ivtx,p)
        entrydata["shower_shlengths"].append( [ showerreco.getVertexShowerShlength(ivtx,p) for p in xrange(3) ] )
        entrydata["vertex_pos"].append( [ showerreco.getVertexPos(ivtx).at(p) for p in xrange(3) ] )
        entrydata["shower_gap"].append( [ showerreco.getVertexShowerGap(ivtx,p) for p in xrange(3) ] )
        entrydata["shower_direction_3d"].append( [ showerreco.getFirstDirection(ivtx,dir) for dir in xrange(3) ])
        entrydata["shower_direction_2d"].append( [ showerreco.getVertexShowerDirection2D(ivtx,dir) for dir in xrange(3) ])
        entrydata["shower_opening_2d"].append( [ showerreco.getVertexShowerOpening2D(ivtx,dir) for dir in xrange(3) ])
        # showerstart also needs a loop over x,y,z
        for p in xrange(3):
            entrydata["shower_start_2d"].append( [ showerreco.getShowerStart2D(ivtx,p,dir) for dir in xrange(3) ])

    # save vertex truth information
    if args.has_mc:
        # pi0 variables...
        entrydata["ccnc"] = []
        entrydata["haspi0"] = []
        entrydata["truefid"] = []
        entrydata["numtrueshowers"] =[]

        for ivtx in xrange(showerreco.numVertices()):
            entrydata["ccnc"].append(showerreco.getCCNC())
            entrydata["haspi0"].append(showerreco.getHasPi0())
            entrydata["truefid"].append(showerreco.getTrueFid())
            entrydata["numtrueshowers"].append(showerreco.getNumTrueShowers())


        if (showerreco.getTrueFid()==1 and (showerreco.getNumTrueShowers() ==1 or showerreco.getNumTrueShowers() ==2)):
            print "Num showers: ",showerreco.getNumTrueShowers()
            entrydata["shower_energy_true"]=[]
            entrydata["shower_recotrue_dist"]=[]
            entrydata["first_direction_true"]=[]
            entrydata["shower_start_2d_true"]=[]

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["shower_energy_true"].append( [ showerreco.getTrueShowerEnergy(ivtx)])
                entrydata["shower_recotrue_dist"].append( [ showerreco.getShowerRecoTrueDist(ivtx)])
                entrydata["first_direction_true"].append( [ showerreco.getTrueShowerDirection(ivtx,dir) for dir in xrange(3)])
                entrydata["shower_start_2d_true"].append( [ showerreco.getTrueShower2DStart(ivtx,idx) for idx in xrange(4)])

        if (showerreco.getTrueFid()==1 and showerreco.getNumTrueShowers()==2):
            entrydata["secondshower_energy_true"]=[]
            entrydata["secondshower_recotrue_dist"]=[]
            entrydata["second_direction_true"]=[]
            entrydata["secondshower_start_2d_true"]=[]

            for ivtx in xrange(showerreco.numVertices()):
                entrydata["secondshower_energy_true"].append( [ showerreco.getTrueSecondShowerEnergy(ivtx)])
                entrydata["secondshower_recotrue_dist"].append( [ showerreco.getSecondShowerRecoTrueDist(ivtx)])
                entrydata["second_direction_true"].append( [ showerreco.getTrueSecondShowerDirection(ivtx,dir) for dir in xrange(3)])
                entrydata["secondshower_start_2d_true"].append( [ showerreco.getTrueSecondShower2DStart(ivtx,idx) for idx in xrange(4)])


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
        # Save second shower output
        entrydata["secondshower_energies"] = []
        entrydata["secondshower_sumQs"] = []
        entrydata["secondshower_shlengths"] = []
        entrydata["secondshower_gap"] = []
        entrydata["pi0mass"] = []
        entrydata["opening_angle_3d"] = []
        entrydata["shower_impact"] = []
        entrydata["secondshower_impact"] =[]
        entrydata["secondshower_direction_3d"]=[]
        entrydata["secondshower_direction_2d"]=[]
        entrydata["secondshower_opening_2d"]=[]
        entrydata["secondshower_start_2d"] =[]

        for ivtx in xrange(showerreco.numVertices()):
            entrydata["secondshower_energies"].append( [ showerreco.getVertexSecondShowerEnergy(ivtx,p) for p in xrange(3) ] )
            entrydata["secondshower_sumQs"].append( [ showerreco.getVertexSecondShowerSumQ(ivtx,p) for p in xrange(3) ] )
            entrydata["secondshower_shlengths"].append( [ showerreco.getVertexSecondShowerShlength(ivtx,p) for p in xrange(3) ] )
            entrydata["secondshower_gap"].append( [ showerreco.getVertexSecondShowerGap(ivtx,p) for p in xrange(3) ] )
            entrydata["pi0mass"].append([showerreco.getPi0Mass(ivtx)])
            entrydata["opening_angle_3d"].append([showerreco.getAlpha(ivtx)])
            entrydata["shower_impact"].append([showerreco.getImpact1(ivtx)])
            entrydata["secondshower_impact"].append([showerreco.getImpact2(ivtx)])
            entrydata["secondshower_direction_2d"].append( [ showerreco.getVertexSecondShowerDirection2D(ivtx,dir) for dir in xrange(3) ])
            entrydata["secondshower_opening_2d"].append( [ showerreco.getVertexSecondShowerOpening2D(ivtx,dir) for dir in xrange(3) ])
            entrydata["secondshower_direction_3d"].append( [ showerreco.getSecondDirection(ivtx,dir) for dir in xrange(3) ])
            # showerstart also needs a loop over x,y,z
            for p in xrange(3):
                entrydata["secondshower_start_2d"].append( [ showerreco.getSecondShowerStart2D(ivtx,p,dir) for dir in xrange(3) ])

    if args.use_bnb:
        entrydata["pi0mass"] = []
        entrydata["haspi0"] = showerreco.getHasPi0()
        entrydata["ccnc"] = showerreco.getCCNC()
        entrydata["disttoint"] = []
        entrydata["impact1"] = []
        entrydata["impact2"] = []
        entrydata["alpha"] = []
        entrydata["firstdirection"] = []
        entrydata["seconddirection"] = []

        mcpg.buildgraph( iolcv, ioll )
        vtx_v = mcpg.findTrackID(-1).start
        entrydata["true_vertex"] = [ vtx_v[i] for i in xrange(3) ] # get the ROOT node
        offset_v = sce.GetPosOffsets( vtx_v[0], vtx_v[1], vtx_v[2] )
        vtx_sce_v = [ vtx_v[0]-offset_v[0]+0.7,
                      vtx_v[1]+offset_v[1],
                      vtx_v[2]+offset_v[2] ]
        entrydata["true_vertex_sce"] = vtx_sce_v

        for ivtx in xrange(showerreco.numVertices()):
                entrydata["pi0mass"].append( showerreco.getPi0Mass(ivtx))
                entrydata["disttoint"].append( showerreco.getDistToInt(ivtx))
                entrydata["impact1"].append( showerreco.getImpact1(ivtx))
                entrydata["impact2"].append( showerreco.getImpact2(ivtx))
                entrydata["alpha"].append( showerreco.getAlpha(ivtx))
                for dir in xrange(3):
                    entrydata["firstdirection"].append( showerreco.getFirstDirection(ivtx,dir))
                    entrydata["seconddirection"].append( showerreco.getSecondDirection(ivtx,dir))


    if args.use_ncpi0:
        entrydata["true_shower_energies"] = []
        entrydata["true_shower_starts"] = []
        entrydata["remaining_adc"] = []
        entrydata["overlap_fraction1"] = []
        entrydata["overlap_fraction2"] = []
        entrydata["purity"] = []
        entrydata["efficiency"] = []
        entrydata["pi0mass"] = []
        entrydata["useformass"] = []
        entrydata["disttoint"] = []
        entrydata["impact1"] = []
        entrydata["impact2"] = []

        for ivtx in xrange(showerreco.numVertices()):
                entrydata["useformass"].append( showerreco.getUseForMass(ivtx))
                entrydata["pi0mass"].append( showerreco.getPi0Mass(ivtx))
                entrydata["disttoint"].append( showerreco.getDistToInt(ivtx))
                entrydata["impact1"].append( showerreco.getImpact1(ivtx))
                entrydata["impact2"].append( showerreco.getImpact2(ivtx))
                entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(ivtx,plane,0) for plane in xrange(2) ] )
                entrydata["overlap_fraction1"].append( [showerreco.getOverlapFraction1(ivtx,plane,1) for plane in xrange(2) ] )
                entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(ivtx,plane,0) for plane in xrange(2) ] )
                entrydata["overlap_fraction2"].append( [showerreco.getOverlapFraction2(ivtx,plane,1) for plane in xrange(2) ] )
                entrydata["purity"].append( [showerreco.getShowerTruthMatchPur(ivtx,shower) for shower in xrange(6)])
                entrydata["efficiency"].append( [showerreco.getShowerTruthMatchEff(ivtx,shower) for shower in xrange(6)])


        for ivtx in xrange(showerreco.numShowers()):
            entrydata["true_shower_energies"].append( [ showerreco.getTrueShowerEnergy(ivtx) for shower in xrange(2) ] )
            entrydata["true_shower_starts"].append( [ showerreco.getTrueShowerStarts(ivtx).at(p) for p in xrange(3) ] )
            entrydata["remaining_adc"].append( [showerreco.getRemainingADC()])

    if args.use_nueint:
        entrydata["uplane_profile"] = []
        entrydata["vplane_profile"] = []
        entrydata["yplane_profile"] = []
        entrydata["purity"] = []
        entrydata["efficiency"] = []


        for ivtx in xrange(showerreco.numVertices()):
            entrydata["purity"].append( [showerreco.getShowerTruthMatchPur(ivtx,shower) for shower in xrange(3)])
            entrydata["efficiency"].append( [showerreco.getShowerTruthMatchEff(ivtx,shower) for shower in xrange(3)])

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

showerreco.finalize()

if args.output_format in ['json','both']:
    fout = open(jout_name, 'w' )
    json.dump( data, fout )
    fout.close()

print "close out"

outll.close()
iolcv.finalize()
if args.input_larlite is not None:
    ioll.close()
