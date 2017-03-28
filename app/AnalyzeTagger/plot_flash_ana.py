import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

rfile = rt.TFile("output_flash_timing.root")

tree = rfile.Get("flashtimes")

nentries = tree.GetEntries()

print "Number of entries in tree: ",nentries

ientry = 0
bytesread = tree.GetEntry(ientry)

plot_types = ["vtxmatched","notvtxmatched","data"]

# DEFINE HISTOGRAMS

# fractional difference from true and hypothesis flash
hpediff = {}
for ptype in plot_types[:2]:
    hpediff[ptype] = rt.TH1D("hpediff_"+ptype,";(PE_{hypo}-PE_{data})/PE_{data}; AU",400,-20.0,20.0)

# hpediff versus pos[0]
hpediff_v_pos = {}
hpediff_v_pos["vtxmatched"] = rt.TH2D("hpediff_v_x_vtxmatched","Vertex-matched ROIs;x (cm);(PE_{hypo}-PE_{data})/PE_{data}",30,-25,275,75,-10.0,20.0)
hpediff_v_pos["notvtxmatched"] = rt.TH2D("hpediff_v_x_notvtxmatched","Not vertex-matched ROIs;x (cm);(PE_{hypo}-PE_{data})/PE_{data}",30,-25,275,75,-10.0,20.0)

# chi-squared
hchi2 = {}
hchi2_v_x = {}
for ptype in plot_types[:2]:
    hchi2[ptype]     = rt.TH1D("hchi2_"+ptype,";#chi^{2}; counts",50,0,1000.0)
    hchi2_v_x[ptype] = rt.TH2D("hchi2_v_x_"+ptype,";x (cm);#chi^{2}",30,-25,275,50,0,1000.0)


while bytesread>0:
    print "entry ",ientry

    # fraction difference between data and hypothesis flash
    for i in range(0,tree.totpe_data.size()):
        totpe_data = tree.totpe_data[i]
        if totpe_data==0:
            continue

        # PE DIFF
        for j in range(0,tree.totpe_hypo_trueflashes.size()):
            pediff = tree.totpe_hypo_trueflashes[j]-totpe_data
            pefracdiff = pediff/totpe_data
            hpediff["vtxmatched"].Fill( pefracdiff )
            hpediff_v_pos["vtxmatched"].Fill( tree.pos[0], pefracdiff )
        for j in range(0,tree.totpe_hypo_falseflashes.size()):
            pediff = tree.totpe_hypo_falseflashes[j]-totpe_data
            pefracdiff = pediff/totpe_data
            hpediff["notvtxmatched"].Fill( pefracdiff )
            hpediff_v_pos["notvtxmatched"].Fill(tree.pos[0],  pefracdiff )

    # CHI2
    for j in range(0,tree.smallest_chi2_trueflashes.size()):
        hchi2["vtxmatched"].Fill( tree.smallest_chi2_trueflashes[j] )
        hchi2_v_x["vtxmatched"].Fill( tree.pos[0], tree.smallest_chi2_trueflashes[j] )
    for j in range(0,tree.smallest_chi2_falseflashes.size()):
        hchi2["notvtxmatched"].Fill( tree.smallest_chi2_falseflashes[j] )
        hchi2_v_x["notvtxmatched"].Fill( tree.pos[0], tree.smallest_chi2_falseflashes[j] )

    ientry += 1
    bytesread = tree.GetEntry(ientry)

# CANVASES

# pediff
cpediff = rt.TCanvas("pediff","",800,600)
hpediff["vtxmatched"].SetLineColor(rt.kRed)
hpediff["notvtxmatched"].Draw()
hpediff["vtxmatched"].Draw("same")
cpediff.Draw()
cpediff.SetLogy(1)
cpediff.Update()

# pediff v pos
cpediff_v_x = rt.TCanvas("pediff_v_pos","PE diff versus x",1400,600)
cpediff_v_x.Divide(2,1)
cpediff_v_x.cd(1)
hpediff_v_pos["vtxmatched"].Draw("COLZ")
cpediff_v_x.cd(2)
hpediff_v_pos["notvtxmatched"].Draw("COLZ")
cpediff_v_x.Draw()
cpediff_v_x.Update()

# Chi2
cchi2 = rt.TCanvas("cchi2","",800,600)
cchi2.Draw()
hchi2["vtxmatched"].SetLineColor(rt.kRed)
hchi2["notvtxmatched"].SetMinimum(0)
hchi2["notvtxmatched"].Draw()
hchi2["vtxmatched"].Draw("same")
cchi2.Update()

# pediff v x
cchi2_v_x = rt.TCanvas("chi2_v_pos","Chi-squared versus x",1400,600)
cchi2_v_x.Divide(2,1)
cchi2_v_x.cd(1)
hchi2_v_x["vtxmatched"].Draw("COLZ")
cchi2_v_x.cd(2)
hchi2_v_x["notvtxmatched"].Draw("COLZ")
cchi2_v_x.Draw()
cchi2_v_x.Update()


raw_input()
