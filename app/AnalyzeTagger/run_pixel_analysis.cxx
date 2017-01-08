#include "Base/DataCoordinator.h"
#include "Base/PSet.h"

// larcv
#include "DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/mctruth.h"

#include "TFile.h"
#include "TTree.h"

int main( int nargs, char** argv ) {

  // run pixel analysis. use 

	// we need to have a data coordinator for each stage because the number of entries could be different.
	// we'll coordinate by using event,subrun,run information
	std::string data_folder = "~/data/larbys/cosmic_tagger/mcc7_bnbcosmic/";

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
	dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv_p1.root", "larcv" ); //stopmu-tagger output


  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_thrumu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_stopmu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();


  // seutp output

  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
  TTree* tree = new TTree("pixelanalysis", "Pixel-level analysis");
  int ncosmic_pixels; // number of non-neutrino pixels
  int nnu_pixels;     // number of neutrino pixels
  int nvertex_pixels; // number of neutrino pixels within some pixel radius of vertex
  int ncosmic_tagged; // number of non-neutrino pixels tagged
  int nnu_tagged;     // number of neutrino pixels tagged
  int nvertex_tagged; // number of neutrino pixels within some pixel radius of vertex tagged
  int mode;           // interaction mode
  int current;        // interaction cufrrent
  float EnuGeV;       // neutrino energy in GeV
  float dwall;        // dwall
  float frac_cosmic;  // fraction of cosmic tagged
  float frac_nu;      // fraction of neutrino pixels tagged
  float frac_vertex;  // fraction of vertex pixels tagged
	tree->Branch("ncosmic_pixels",&ncosmic_pixels,"ncosmic_pixels/I");
	tree->Branch("nnu_pixels",&nnu_pixels,"nnu_pixels/I");
	tree->Branch("nvertex_pixels",&nvertex_pixels,"nvertex_pixels/I");
  

  int nentries = dataco_stopmu.get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

  	dataco_stopmu.goto_entry(ientry,"larcv");

  	int run, subrun, event;
  	dataco_stopmu.get_id(run,subrun,event);

  	dataco_thrumu.goto_event(run,subrun,event,"larcv");
  	dataco_source.goto_event(run,subrun,event,"larcv");

  	// ok now to do damage


  	break;

  }


  return 0;
}
