import os,sys
import ROOT as rt

if len(sys.argv)!=3:
    print "use: plot_example.py [dim] [val]"

rt.gStyle.SetOptStat(0)
rt.gStyle.SetPadRightMargin(0.2)

rfile = rt.TFile("reverse_sce_table.root", "OPEN")
xorigin = rfile.Get("xorigin")
yorigin = rfile.Get("yorigin")
zorigin = rfile.Get("zorigin")
horigin = [ xorigin, yorigin, zorigin ]

dim = str(sys.argv[1].strip())
value = float(sys.argv[2].strip())

print "dim=",dim,"value=",value

haxis = None
hists = {}
if "x" in dim:
    haxis = xorigin.GetXaxis()
    hxbin = haxis.FindBin( value )
    ybins = xorigin.GetYaxis().GetNbins()
    ymin = xorigin.GetYaxis().GetXmin()
    ymax = xorigin.GetYaxis().GetXmax()
    zbins = xorigin.GetZaxis().GetNbins()
    zmin = xorigin.GetZaxis().GetXmin()
    zmax = xorigin.GetZaxis().GetXmax()
    hists[0] = rt.TH2D( "hyz_xorigin", "x-origin", ybins, ymin, ymax, zbins, zmin, zmax )
    hists[1] = rt.TH2D( "hyz_yorigin", "y-origin", ybins, ymin, ymax, zbins, zmin, zmax )
    hists[2] = rt.TH2D( "hyz_zorigin", "z-origin", ybins, ymin, ymax, zbins, zmin, zmax )

    for i in range(0,3):
        for ybin in range(1,ybins+1):
            for zbin in range(1,zbins+1):
                val = horigin[i].GetBinContent( hxbin, ybin, zbin )
                hists[i].SetBinContent( ybin, zbin, val )
            
elif dim=="y":
    haxis = xorigin.GetYaxis()
elif dim=="z":
    haxis = xorigin.GetZaxis()

c = rt.TCanvas("c","c",1400,500)
c.Draw()
c.Divide(3,1)
for i in range(0,3):
    c.cd(i+1)
    hists[i].Draw("COLZ")

c.Update()

print "done"
raw_input()
