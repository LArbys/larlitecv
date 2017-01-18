import os,sys
import ROOT as rt
from math import sqrt

rfilename = sys.argv[1]
rfile = rt.TFile(rfilename)

tree = rfile.Get("flashroi")

nentries = tree.GetEntries()
chi2_cut_value = 50.0
scale_cut_value = 0.1
stop_each_event = True

c = rt.TCanvas("pespectrum","PE Spectrum", 12000,600)
c.Divide(2,1)
c.Draw()
hmeasured  = rt.TH1D("hmeasured","Measured",32,0,32)
hpredicted = rt.TH1D("hpredicted","Predicted", 32, 0, 32)
hinit      = rt.TH1D("hinit","Initial Prediction", 32, 0, 32)

hmeasured.SetLineWidth(3)
hpredicted.SetLineColor(rt.kRed)
hinit.SetLineColor(rt.kGreen+3)

hscale_nu = rt.TH1D("hscale_nu","",200,0,10)
hscale_bg = rt.TH1D("hscale_bg","",200,0,10)
hscale_uf = rt.TH1D("hscale_uf","",200,0,10)
hscale_bg.SetLineColor(rt.kRed)
hscale_uf.SetLineColor(rt.kGreen+3)

hchi2_fit_nu = rt.TH1D("hchi2_fit_nu","",200,0,100)
hchi2_fit_bg = rt.TH1D("hchi2_fit_nu","",200,0,100)
hchi2_raw_nu = rt.TH1D("hchi2_raw_nu","",200,0,100)
hchi2_raw_bg = rt.TH1D("hchi2_raw_bg","",200,0,100)
hchi2_fit_nu.SetLineWidth(2)
hchi2_fit_nu.SetLineColor(rt.kRed)
hchi2_fit_bg.SetLineColor(rt.kRed)
hchi2_raw_nu.SetLineWidth(2)
hchi2_raw_nu.SetLineColor(rt.kGreen+4)
hchi2_raw_bg.SetLineColor(rt.kGreen+4)

h2d = rt.TH2D("h2d","",50,0,200,50,0,1.0)
hchi2_v_scale_nu = rt.TH2D("hchi2_v_scale_nu",";Chi2 vs. Scale", 100,0,10.0, 50,0,100 )
hchi2_v_scale_bg = rt.TH2D("hchi2_v_scale_bg",";Chi2 vs. Scale", 100,0,10.0, 50,0,100 )
hroi_per_event = rt.TH1D("hroi_per_event","",50,0,50)

record = open('output_chi2.txt','w')

npass_nu = 0
npass_bg = 0

nrejected_nu = 0
nrejected_bg = 0

nevent_rois = 0
current_rse = (0,0,0)
predicted_hists = {}
init_hists = {}
prediction_idx = 0
for i in range(0,20):
	predicted_hists[i] = hpredicted.Clone("hprediction_%d"%(i))
	init_hists[i] = hpredicted.Clone("hinit_%d"%(i))

for ientry in range(0,nentries):
	tree.GetEntry(ientry)

	entry_rse = ( tree.run, tree.subrun, tree.event )
	if entry_rse!=current_rse:
		if current_rse!=(0,0,0):
			print "about to start new event"
			c.Update()
                        if stop_each_event:
			        raw_input()
                hroi_per_event.Fill( nevent_rois )
                nevent_rois = 0
		print "========================================"
		print "New Event: [%d, %d, %d]" % ( entry_rse )
		hmeasured.Reset()
		for idx,hist in predicted_hists.items():
			hist.Reset()
			hist.SetLineColor(rt.kRed)
			hist.SetLineWidth(1)
			hist.SetLineStyle(2)
		for idx,hist in init_hists.items():
			hist.Reset()
			hist.SetLineColor(rt.kGreen+3)
			hist.SetLineWidth(1)	
		current_rse = entry_rse		
		prediction_idx = 0
		filled_measured = False
		c.cd(1).Clear()
		c.cd(2).Clear()

	if prediction_idx>=20:
		continue

	print "  [entry %d]"%(ientry)
	print "  flash range (wires): [ %.2f, %.2f ]" % ( tree.flash_range[0], tree.flash_range[1] )
	print "  flash range (z):     [ %.2f, %.2f ]" % ( tree.flash_range[0]*0.3, tree.flash_range[1]*0.3 )	
	print "  nupixelfrac=%.3f"%(tree.fracnupixels)
	print "  score: ",tree.score

	chi2_fit = 0.;
	chi2_raw = 0.;
	raw_total = tree.rawtotpe;
	for ipmt in range(0,32):
		if not filled_measured:
			hmeasured.SetBinContent(ipmt+1, tree.measured[ipmt])
		else:
			filled_measured = True
		predicted_hists[prediction_idx].SetBinContent(ipmt+1, tree.hypothesis[ipmt])
		init_hists[prediction_idx].SetBinContent(ipmt+1, tree.raw_hypothesis[ipmt])

		if tree.measured[ipmt]+tree.hypothesis[ipmt]>0:
			diff = tree.measured[ipmt]-tree.hypothesis[ipmt]
			err = sqrt( tree.measured[ipmt] + tree.hypothesis[ipmt] )
			chi2_fit += diff*diff/(err*err)
		if tree.measured[ipmt]+tree.raw_hypothesis[ipmt]>0:
			diff = tree.measured[ipmt]-tree.raw_hypothesis[ipmt]
			err = sqrt( tree.measured[ipmt] + tree.raw_hypothesis[ipmt] )
			chi2_raw += diff*diff/(err*err)
	print "  measured pe=%.2f vs. prediceted pe=%.2f : ratio=%.2f" % ( tree.totalpe, tree.predictpe, float(tree.predictpe)/float(tree.totalpe) )
	print "  measured pe=%.2f vs. raw predicted pe=%.2f : ratio=%.2f" % ( tree.totalpe, raw_total, float(raw_total)/float(tree.totalpe) )

	chi2_fit /= 32.0
	chi2_raw /= 32.0
        raw_scale = float(raw_total)/float(tree.totalpe)
	print "  chi2 fit=",chi2_fit," raw=",chi2_raw
	if tree.fracnupixels>0.3:
		predicted_hists[prediction_idx].SetLineWidth(2)
		init_hists[prediction_idx].SetLineWidth(2)
		predicted_hists[prediction_idx].SetLineStyle(1)
		init_hists[prediction_idx].SetLineStyle(1)
		hchi2_fit_nu.Fill( chi2_fit )
		hchi2_raw_nu.Fill( chi2_raw )
	        hchi2_v_scale_nu.Fill( raw_scale, chi2_raw )                
		if chi2_raw<chi2_cut_value and raw_scale>scale_cut_value:
			npass_nu += 1
                        nevent_rois+=1                        
		else:
			nrejected_nu += 1
	else:		
		hchi2_fit_bg.Fill( chi2_fit )
		hchi2_raw_bg.Fill( chi2_raw )
	        hchi2_v_scale_bg.Fill( raw_scale, chi2_raw )
		if chi2_raw<chi2_cut_value and raw_scale>scale_cut_value:                
			npass_bg += 1
                        nevent_rois+=1                        
		else:
			nrejected_bg += 1
	h2d.Fill( chi2_raw, tree.fracnupixels )

	prediction_idx += 1


	c.cd(1)
	hmeasured.Draw()
	for i in range(0,prediction_idx):
		predicted_hists[i].Draw("same")
		init_hists[i].Draw("same")

	if tree.fracnupixels>0.3:
		hscale_nu.Fill( tree.predictpe/tree.totalpe )
		hscale_uf.Fill( init_hists[prediction_idx-1].Integral()/tree.totalpe )
		print >> record,(tree.run,tree.subrun,tree.event)," ",tree.fracnupixels," ",chi2_raw," ",raw_scale
	else:
		hscale_bg.Fill( tree.predictpe/tree.totalpe )

	c.cd(2)
	if hchi2_fit_bg.GetMaximum()>hchi2_fit_nu.GetMaximum():
		hchi2_fit_bg.SetMinimum(0)
		hchi2_fit_bg.Draw()
		hchi2_fit_nu.Draw("same")
		hchi2_raw_bg.Draw("same")
		hchi2_raw_nu.Draw("same")
	else:
		hchi2_raw_bg.SetMinimum(0)		
		hchi2_raw_bg.Draw()
		hchi2_raw_nu.Draw("same")
		hchi2_fit_bg.Draw("same")
		hchi2_fit_nu.Draw("same")


	#if hscale_nu.GetMaximum()>hscale_bg.GetMaximum():
	#	hscale_nu.Draw()
	#	hscale_bg.Draw("same")
	#	hscale_uf.Draw("same")
	#else:
	#	hscale_bg.Draw()
	#	hscale_nu.Draw("same")		
	#	hscale_uf.Draw("same")		

	c.Update()


print "Done with event loop"
record.close()

print "NUM NU: ",npass_nu+nrejected_nu
print "NU PASSED: ",float(npass_nu)/float(npass_nu+nrejected_nu)
print "BG REJECTED: ",float(nrejected_bg)/float(npass_bg+nrejected_bg)

c2 = rt.TCanvas("c2","c2",1200,600)
c2.Divide(2,1)
c2.cd(1)
c2.Draw()
#h2d.Draw("COLZ")
hchi2_v_scale_bg.Draw("COLZ")
c2.cd(2)
hchi2_v_scale_nu.Draw("COLZ")

c2.Update()

c3 = rt.TCanvas("c3","c3",800,600)
c3.Draw()
hroi_per_event.Draw()
c3.Update()

raw_input()



