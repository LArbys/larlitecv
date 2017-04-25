import os,sys
import ROOT as rt

if len(sys.argv)!=2:
    print "usage: python plot_tagger_ana.py [pixel analysis root file]"
    print "this root file is made using run_pixel_analysis"
    sys.exit(-1)

rt.gStyle.SetOptStat(0)

#rfile = rt.TFile("output_pixel_analysis.root")
#rfile = rt.TFile("output_pixel_analysis_rerunpad20.root")
#rfile = rt.TFile("pixana_extbnb.root")
rfile = rt.TFile(sys.argv[1])

tree = rfile.Get("pixana")

nentries = tree.GetEntries()


print "Number of entries in tree: ",nentries

ientry = 0
bytesread = tree.GetEntry(ientry)

stages = ["thrumu","stopmu","untagged","croi"]
colors = {"thrumu":rt.kBlack,
          "stopmu":rt.kBlue,
          "untagged":rt.kRed,
          "croi":rt.kGreen-3}

# DEFINE HISTOGRAMS

# fraction cosmic+vertex pixels tagged
htagged = {}
for stage in stages:
    htagged[stage]= rt.TH1D("htagged_"+stage,";fraction of pixels; counts",50,0,1.01)
    htagged[stage].SetLineColor( colors[stage] )
    htagged[stage].SetLineWidth(2)

# fraction cosmic+vertex tagged versus total pixels (complexity)
htagged_v_totpixels = {}
for stage in stages:
    htagged_v_totpixels[stage]= rt.TH2D("htagged_"+stage,stage+";total pixels (1k); fraction tagged",50,0,1e5,50,0,1.01)
    htagged_v_totpixels[stage].SetLineColor( colors[stage] )
    htagged_v_totpixels[stage].SetLineWidth(2)

# fraction cosmic+vertex pixels tagged
hvertex = {}
for stage in stages:
    hvertex[stage]= rt.TH1D("hvertex_"+stage,";fraction of vertex pixels; counts",50,0,1.01)
    hvertex[stage].SetLineColor( colors[stage] )
    hvertex[stage].SetLineWidth(2)

# fraction cosmic+vertex pixels tagged
hvertex_badch = rt.TH1D("hvertex_badch",";fraction of vertex pixels are bad",50,0,1.01)
hvertex_badch_incroi = rt.TH1D("hvertex_badch_incroi",";fraction of vertex pixels are bad",50,0,1.01)
hvertex_badch_nocroi = rt.TH1D("hvertex_badch_nocroi",";fraction of vertex pixels are bad",50,0,1.01)
hvertex_badch.SetLineWidth(2)
    

# fraction vertex tagged versus total pixels (complexity)
hvertex_v_totpixels = {}
for stage in stages:
    hvertex_v_totpixels[stage]= rt.TH2D("hvertex_v_totpixels"+stage,stage+";total pixels (1k); fraction of vertex pixels tagged",50,0,1e5,50,0,1.01)
    hvertex_v_totpixels[stage].SetLineColor( colors[stage] )
    hvertex_v_totpixels[stage].SetLineWidth(2)

# fraction vertex tagged versus total pixels (complexity)
hvertex_v_dwall = {}
for stage in stages:
    hvertex_v_dwall[stage]= rt.TH2D("hvertex_v_dwall"+stage,stage+";distance to wall (cm); fraction of vertex pixels tagged",50,0,100,50,0,1.001)
    hvertex_v_dwall[stage].SetLineColor( colors[stage] )
    hvertex_v_dwall[stage].SetLineWidth(2)

# fraction cosmic+vertex tagged only one versus total pixels (complexity)
htag1_v_totpixels = {}
for stage in stages:
    htag1_v_totpixels[stage]= rt.TH2D("htag1_"+stage,stage+";total pixels (1k); fraction tagged only once",50,0,1e5,50,0,1.0)
    htag1_v_totpixels[stage].SetLineColor( colors[stage] )
    htag1_v_totpixels[stage].SetLineWidth(2)

# fraction cosmic+vertex tagged many times versus total pixels (complexity)
htag2_v_totpixels = {}
for stage in stages:
    htag2_v_totpixels[stage]= rt.TH2D("htag2_"+stage,stage+";total pixels (1k); fraction tagged more than once",50,0,1e5,50,0,1.0)
    htag2_v_totpixels[stage].SetLineColor( colors[stage] )
    htag2_v_totpixels[stage].SetLineWidth(2)

# num rois
hnrois = rt.TH1D("hnroi",";number of ROIs per event",50,0,50)

# fraction of true crossing points tagged
hendtagged = rt.TH1D("hendtagged",";fraction of true end points tagged per event",50,0,1.0001)
hendtagged_v_totpix = rt.TH2D("hendtagged",";total pixels;fraction of end points tagged",50,0,1e5,50,0,1.0)

# efficiency versus position
hpos = {}
hpos_vxt_incroi = {}
hratio = {}
axes = {0:"x",1:"y",2:"z"}
hpos[0] = rt.TH1D("hpos_x_all","; x vertex pos (cm)",15,-25,275)
hpos[1] = rt.TH1D("hpos_y_all","; x vertex pos (cm)",20,-120,120)
hpos[2] = rt.TH1D("hpos_z_all","; x vertex pos (cm)",20,-50,1050)
for i in range(0,3):
    hpos_vxt_incroi[i] = hpos[i].Clone( str(hpos[i].GetName()).replace("all","incroi") )
    hratio[i] = hpos[i].Clone( str(hpos[i].GetName()).replace("all","ratio") )

# efficiency versus total pixels
hefftot_all    = rt.TH1D("hefftot_all",";total pixels; efficiency",10,0,1e5)
hefftot_incroi = hefftot_all.Clone( str(hefftot_all.GetName()).replace("all","incroi") )
hefftot_ratio  = hefftot_all.Clone( str(hefftot_all.GetName()).replace("all","ratio") )

while bytesread>0:
    print "entry ",ientry

    # cosmic tagging
    tottagged = 0
    tagged_many = 0
    tagged_once = 0
    vertex_tagged = 0
    totpixels = 0
    for idx,stage in enumerate(stages[:3]):
        tottagged   += tree.ncosmic_tagged[4*idx+3] + tree.nvertex_tagged[4*idx+3]
        tagged_many += tree.ncosmic_tagged_many[4*idx+3]
        tagged_once += tree.ncosmic_tagged_once[4*idx+3]
        vertex_tagged += tree.nvertex_tagged[4*idx+3]
        totpixels = tree.nthreshold_pixels[3]
        if totpixels>0:
            frac = float(tottagged)/float(totpixels)
            frac_many = float(tagged_many)/float(totpixels)
            frac_once = float(tagged_once)/float(totpixels)
        else:
            frac = 0.0
            frac_many = 0.0
            frac_once = 0.0
        if tree.nvertex_pixels[0]>0:
            vertex_frac = float(vertex_tagged)/float(tree.nvertex_pixels[3])
        else:
            vertex_frac = 0.0
        htagged[stage].Fill( frac )
        hvertex[stage].Fill( vertex_frac )
        htagged_v_totpixels[stage].Fill(totpixels,frac)
        hvertex_v_totpixels[stage].Fill(totpixels,vertex_frac)
        hvertex_v_dwall[stage].Fill(tree.dwall,vertex_frac)
        htag1_v_totpixels[stage].Fill(totpixels,frac_once)
        htag2_v_totpixels[stage].Fill(totpixels,frac_many)

    # neutrino tagging
    if totpixels>0 and tree.nvertex_pixels[3]>0:
        frac_vertex = float(tree.nvertex_tagged[4*3+3])/float(tree.nvertex_pixels[3])
        htagged["croi"].Fill( frac_vertex )
        htagged_v_totpixels["croi"].Fill(totpixels,frac_vertex)
        hvertex_v_totpixels["croi"].Fill(totpixels,frac_vertex)


    for i in range(0,3):
        if tree.vtx_in_croi==1:
            hpos_vxt_incroi[i].Fill( tree.pos[i] )
            hratio[i].Fill( tree.pos[i] )
        hpos[i].Fill(tree.pos[i])

    if tree.vtx_in_croi==1:
        hefftot_incroi.Fill( totpixels )
        hefftot_ratio.Fill( totpixels )
    hefftot_all.Fill( totpixels )

    # number of ROIs
    hnrois.Fill( tree.num_rois )

    # fraction of end points tagged
    if tree.true_crossingpoints>0:
        frac_endpoints = float(tree.tagged_crossingpoints)/float(tree.true_crossingpoints)
    else:
        frac_endpoints = 0.
    hendtagged.Fill( frac_endpoints )
    hendtagged_v_totpix.Fill( totpixels, frac_endpoints )

    # fraction of vertex pixels are bad
    if tree.nvertex_pixels[3]>0:
        fracbad = float(tree.nvertex_badch[3])/float(tree.nvertex_pixels[3])
        hvertex_badch.Fill( fracbad )
        if tree.vtx_in_croi==1:
            hvertex_badch_incroi.Fill( fracbad )
        else:
            hvertex_badch_nocroi.Fill( fracbad )

    ientry += 1
    bytesread = tree.GetEntry(ientry)

# CANVASES

# tagged
ctagged = rt.TCanvas("ctagged","Fraction of pixels tagged",1400,600)
#ctagged.Divide(2,1)
ctagged.Draw()
#ctagged.cd(1)
htagged["untagged"].Draw()
htagged["stopmu"].Draw("same")
htagged["thrumu"].Draw("same")
#ctagged.cd(2)
#htagged["croi"].Draw()

# vertex
cvertex = rt.TCanvas("cvertex","Fraction of pixels vertex",1400,600)
#cvertex.Divide(2,1)
cvertex.Draw()
#cvertex.cd(1)
hvertex["untagged"].Draw()
hvertex["stopmu"].Draw("same")
hvertex["thrumu"].Draw("same")
#cvertex.cd(2)
#hvertex["croi"].Draw()

# cosmic tagged versus complexity
ctagged_v_totpixels = rt.TCanvas("ctagged_v_totpixels","Fraction of pixels tagged versus total pixels (complexity)",1200,1200)
ctagged_v_totpixels.Divide(2,2)
ctagged_v_totpixels.Draw()
ctagged_v_totpixels.cd(1)
htagged_v_totpixels["thrumu"].Draw("COLZ")
ctagged_v_totpixels.cd(2)
htagged_v_totpixels["stopmu"].Draw("COLZ")
ctagged_v_totpixels.cd(3)
htagged_v_totpixels["untagged"].Draw("COLZ")

# vertex versus complexity
cvertex_v_totpixels = rt.TCanvas("cvertex_v_totpixels","Fraction of vertex pixels tagged versus total pixels (complexity)",1200,1200)
cvertex_v_totpixels.Divide(2,2)
cvertex_v_totpixels.Draw()
cvertex_v_totpixels.cd(1)
hvertex_v_totpixels["thrumu"].Draw("COLZ")
cvertex_v_totpixels.cd(2)
hvertex_v_totpixels["stopmu"].Draw("COLZ")
cvertex_v_totpixels.cd(3)
hvertex_v_totpixels["untagged"].Draw("COLZ")
cvertex_v_totpixels.cd(4)
hvertex_v_totpixels["croi"].Draw("COLZ")

# vertex versus complexity
cvertex_v_dwall = rt.TCanvas("cvertex_v_dwall","Fraction of vertex pixels tagged versus dwall",1200,1200)
cvertex_v_dwall.Divide(2,2)
cvertex_v_dwall.Draw()
cvertex_v_dwall.cd(1)
hvertex_v_dwall["thrumu"].Draw("COLZ")
cvertex_v_dwall.cd(2)
hvertex_v_dwall["stopmu"].Draw("COLZ")
cvertex_v_dwall.cd(3)
hvertex_v_dwall["untagged"].Draw("COLZ")
cvertex_v_dwall.cd(4)
hvertex_v_dwall["croi"].Draw("COLZ")

# many tagged versus complexity
ctag1_v_totpixels = rt.TCanvas("ctag1_v_totpixels","Fraction of pixels tagged many times versus total pixels (complexity)",1200,1200)
ctag1_v_totpixels.Divide(2,2)
ctag1_v_totpixels.Draw()
ctag1_v_totpixels.cd(1)
htag1_v_totpixels["thrumu"].Draw("COLZ")
ctag1_v_totpixels.cd(2)
htag1_v_totpixels["stopmu"].Draw("COLZ")
ctag1_v_totpixels.cd(3)
htag1_v_totpixels["untagged"].Draw("COLZ")

# many tagged versus complexity
ctag2_v_totpixels = rt.TCanvas("ctag2_v_totpixels","Fraction of pixels tagged many times versus total pixels (complexity)",1200,1200)
ctag2_v_totpixels.Divide(2,2)
ctag2_v_totpixels.Draw()
ctag2_v_totpixels.cd(1)
htag2_v_totpixels["thrumu"].Draw("COLZ")
ctag2_v_totpixels.cd(2)
htag2_v_totpixels["stopmu"].Draw("COLZ")
ctag2_v_totpixels.cd(3)
htag2_v_totpixels["untagged"].Draw("COLZ")

# number of rois per event
rt.gStyle.SetOptStat(1)
cnrois = rt.TCanvas("cnrois","Number of ROIs",800,600)
cnrois.Draw()
hnrois.SetLineWidth(2)
hnrois.Draw()
rt.gStyle.SetOptStat(0)

# vertex distributions
cvtxpos = rt.TCanvas("cvtxpos","True Vertex Distributions",800,1600)
cvtxpos.Divide(2,3)
for i in range(0,3):
    cvtxpos.cd( 2*i+1 )
    hpos[i].Draw()
    hpos_vxt_incroi[i].SetLineColor(rt.kRed)
    hpos_vxt_incroi[i].Draw("same")
    hratio[i].Divide( hpos[i] )
    cvtxpos.cd(2*i+2)
    hratio[i].Draw()

# efficiency versus total number of pixels
cefftotpixels = rt.TCanvas("cefftotpixels","Efficiency versus total number of pixels",1400,600)
cefftotpixels.Divide(2,1)
cefftotpixels.cd(1)
hefftot_all.Draw()
hefftot_incroi.SetLineColor(rt.kRed)
hefftot_incroi.Draw("same")
cefftotpixels.cd(2)
hefftot_ratio.Divide( hefftot_all )
hefftot_ratio.Draw()

# number of end points tagged
cendpoints = rt.TCanvas("cendpoints","Number of End Points Tagged", 1400,600)
cendpoints.Divide(2,1)
cendpoints.cd(1)
hendtagged.SetLineWidth(2)
hendtagged.Draw()
cendpoints.cd(2)
hendtagged_v_totpix.Draw("COLZ")

# fraction bad channels
cfracbad = rt.TCanvas("cbadch","Fraction of vertex pixels bad",1400,600)
cfracbad.Divide(2,1)
cdf_vertex_badch = hvertex_badch.Clone("cdf_vertex_badch")
cdf_vertex_badch.Reset()
cdf_vertex_badch.SetTitle(";fraction of vertex pixels are bad;fraction of events")
tot = 0.
if hvertex_badch.Integral()>0:
    for b in range(1,hvertex_badch.GetXaxis().GetNbins()+1):
        tot += hvertex_badch.GetBinContent(b)
        cdf_vertex_badch.SetBinContent(b,tot/hvertex_badch.Integral())
cfracbad.cd(1)
hvertex_badch_incroi.SetLineColor(rt.kRed)
hvertex_badch_nocroi.SetLineColor(rt.kBlack)
hvertex_badch.Draw()
hvertex_badch_incroi.Draw("same")
hvertex_badch_nocroi.Draw("same")
cfracbad.cd(2)
cdf_vertex_badch.Draw()
cfracbad.Draw()

print "Mean: ",hnrois.GetMean()
print "Efficiency: ",hpos_vxt_incroi[0].Integral()/hpos[0].Integral()

raw_input()
