import os,sys

input_larlite = sys.argv[1]

import ROOT as rt
from larlite import larlite
from larlitecv import larlitecv

io = larlite.storage_manager( larlite.storage_manager.kREAD )
io.add_in_filename( input_larlite )
io.open()

outfile = rt.TFile("out_larbysmc.root","recreate");
outtree = rt.TTree("larbysmc", "Exracted truth quantities for DL LEE analyses")

lmc = larlitecv.LArbysMC()

outfile.cd()
lmc.bindAnaVariables(outtree)

nentries = io.get_entries()
for ientry in xrange(nentries):
    print "[[ ENTRY ",ientry," ]]"
    io.go_to(ientry)
    lmc.process(io)
    lmc.printInteractionInfo()
    outtree.Fill()

outfile.Write()
io.close()
print "[FIN]"
