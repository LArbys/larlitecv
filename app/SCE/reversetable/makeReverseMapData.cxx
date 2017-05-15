#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH3D.h"

#include "SpaceChargeMicroBooNE.h"

int main () {
  std::cout << "Make Reverse Map Data" << std::endl;

  larlitecv::SpaceChargeMicroBooNE sce;
  double gridranges[3][2] = { {-10.0, 260.0},
			      {-120.0,120.0},
			      {0.0, 1050.0} };
  int nbins[3]  = { 54, 48, 210 };
  //int nbins[3]  = { 27, 24, 105 };  
  int nsteps[3] = { 0 };
  float stepsize[3] = {0};
  for (int i=0; i<3; i++) {
    nsteps[i] = nbins[i]*3;
    stepsize[i] = (gridranges[i][1]-gridranges[i][0])/float(nsteps[i]);
  }

  int totsteps = (nsteps[0]+1)*(nsteps[1]+1)*(nsteps[2]+1);

  int iter = 0;

  TFile* rfile = new TFile( "reverse_sce_table.root", "RECREATE" );
  TH3D hentries("hentries","",
		nbins[0], gridranges[0][0], gridranges[0][1],
		nbins[1], gridranges[1][0], gridranges[1][1],
		nbins[2], gridranges[2][0], gridranges[2][1] );
  TH3D hxorigin("xorigin","",
		nbins[0], gridranges[0][0], gridranges[0][1],
		nbins[1], gridranges[1][0], gridranges[1][1],
		nbins[2], gridranges[2][0], gridranges[2][1] );
  TH3D hyorigin("yorigin","",
		nbins[0], gridranges[0][0], gridranges[0][1],
		nbins[1], gridranges[1][0], gridranges[1][1],
		nbins[2], gridranges[2][0], gridranges[2][1] );
  TH3D hzorigin("zorigin","",
		nbins[0], gridranges[0][0], gridranges[0][1],
		nbins[1], gridranges[1][0], gridranges[1][1],
		nbins[2], gridranges[2][0], gridranges[2][1] );


  unsigned long totiters = totsteps;
  for (int ix=0; ix<nsteps[0]; ix++) {
    for (int iy=0; iy<nsteps[1]; iy++) {
      for (int iz=0; iz<nsteps[2]; iz++) {

	float x = gridranges[0][0] + (float(ix)+0.5)*stepsize[0];
	float y = gridranges[1][0] + (float(iy)+0.5)*stepsize[1];
	float z = gridranges[2][0] + (float(iz)+0.5)*stepsize[2];

	std::vector<double> offsets = sce.GetPosOffsets( x, y, z);

	float x1 = (x - offsets[0] + 0.17);
	float y1 = (y + offsets[1]);
	float z1 = (z + offsets[2]);


	float dx = x - x1;
	float dy = y - y1;
	float dz = z - z1;
	//float dx = x;
	//float dy = y;
	//float dz = z;

	int xbin = hxorigin.GetXaxis()->FindBin( x1 );
	int ybin = hyorigin.GetYaxis()->FindBin( y1 );
	int zbin = hzorigin.GetZaxis()->FindBin( z1 );	
	
	hentries.Fill( x1, y1, z1 );
	hxorigin.SetBinContent( xbin, ybin, zbin, hxorigin.GetBinContent(xbin, ybin, zbin)+dx );
	hyorigin.SetBinContent( xbin, ybin, zbin, hyorigin.GetBinContent(xbin, ybin, zbin)+dy );
	hzorigin.SetBinContent( xbin, ybin, zbin, hzorigin.GetBinContent(xbin, ybin, zbin)+dz );
	
	iter++;

	if ( iter%100000==0 ) {
	  std::cout << "iter " << iter << " of " << totiters << " " << float(iter)/float(totiters)*100.0 << "% "
		    << " dx=" << dx
		    << "pos=(" << x << "," << y << "," << z << ") "
		    << "pos\'=(" << x1 << "," << y1 << "," << z1 << ") "
		    << "bins=(" << xbin << "," << ybin << "," << zbin << ") "
		    <<  std::endl;
	}
      }
    }	
  }

  std::cout << "Divide" << std::endl;
  int nfilled = 0;
  int totbins = 0;
  for (int xbin=1; xbin<=hentries.GetXaxis()->GetNbins(); xbin++) {
    for (int ybin=1; ybin<=hentries.GetYaxis()->GetNbins(); ybin++) {
      for (int zbin=1; zbin<=hentries.GetZaxis()->GetNbins(); zbin++) {
	totbins++;
	int nentries = hentries.GetBinContent( xbin, ybin, zbin );
	if (nentries==0)
	  continue;
	hxorigin.SetBinContent( xbin, ybin, zbin, hxorigin.GetBinContent( xbin, ybin, zbin )/float(nentries) );
	hyorigin.SetBinContent( xbin, ybin, zbin, hyorigin.GetBinContent( xbin, ybin, zbin )/float(nentries) );
	hzorigin.SetBinContent( xbin, ybin, zbin, hzorigin.GetBinContent( xbin, ybin, zbin )/float(nentries) );
	nfilled++;
      }
    }
  }
  std::cout << "fraction filled: " << float(nfilled)/float(totbins) << std::endl;

  std::cout << "Save" << std::endl;
  //hxorigin.Write();
  //hyorigin.Write();
  //hzorigin.Write();
  rfile->Write();
  
  return 0;
}
