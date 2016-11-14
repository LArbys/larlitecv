import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

fgood = open("spoon_roi_scan.csv",'r')
lgood = fgood.readlines()
use_handscan = False
use_dwallcut = False

goodroi = []
for l in lgood[1:]:
    u = int(l.split(",")[1].strip())
    v = int(l.split(",")[2].strip())
    y = int(l.split(",")[3].strip())
    goodroi.append( [u,v,y] )

fin = rt.TFile("output_analysis.root","open")
bmt = fin.Get("bmt")

hnu  = rt.TH1D("hnu","",50,0,1.001)
htot = rt.TH1D("htot","",50,0,1.001)

nentries = bmt.GetEntries()
nbelow90 = 0
for i in range(0,nentries):
    bmt.GetEntry(i)
    if use_dwallcut and bmt.dwall<5:
        continue
    for p in range(0,3):
        totfrac = bmt.total_frac_remain[p]
        nufrac   = bmt.nubb_frac_remain[p]
        if not use_handscan or goodroi[i][p]==1:
            if False or ( bmt.enu<500 and bmt.enu>100 and bmt.dwall>5):
                htot.Fill( totfrac )
                hnu.Fill( nufrac )
                if nufrac<0.9:
                    print i,(bmt.run,bmt.subrun,bmt.event),p,nufrac
                    nbelow90+=1
print "below 90: ",float(nbelow90)/float(hnu.Integral())

c = rt.TCanvas("c","c",800,400)
c.Draw()
hnu.SetLineColor(rt.kRed)
hnu.Draw()
htot.Draw("same")

c.Update()

print "[enter] to exit"
raw_input()
