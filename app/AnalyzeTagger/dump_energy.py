import ROOT as rt

rfile = rt.TFile("output_pixel_analysis.root")
ttree = rfile.Get("pixana")

print ttree

fout = open("event_info.txt",'w')

print >> fout, "run",'\t',"subrun",'\t',"event","\t","dwall vtx",'\t',"energy",'\t',"dwall_lepton",'\t',"lepton boundary",'\t',"num protons >60 MeV",'\t',"prim proton ke"
nentries = ttree.GetEntries();

for ientry in range(0,nentries):
    ttree.GetEntry(ientry)

    print >> fout, ttree.run,"\t",ttree.subrun,"\t",ttree.event,"\t",ttree.dwall,"\t",ttree.EnuGeV,'\t',
    print >> fout, ttree.dwall_lepton,'\t',ttree.lepton_boundary,'\t',ttree.num_protons_over60mev,'\t',ttree.primary_proton_ke


fout.close()
