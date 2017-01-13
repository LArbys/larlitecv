#include <iostream>
#include <string>

#include "Base/DataCoordinator.h"

#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

int main( int nargs, char** argv ) {

	std::cout << "[[ ROI SELECTION ]]" << std::endl;

	// we need to have a data coordinator for each stage because the number of entries could be different.
	// we'll coordinate by using event,subrun,run information
	std::string data_folder = "~/working/data/larbys/cosmic_tagger_dev/";

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
	dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv_p1.root", "larcv" ); //stopmu-tagger output


  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_thrumu.configure( "config.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_stopmu.configure( "config.cfg", "StorageManager", "IOManager", "ContainedROI" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet pixana_cfg = cfg.get<larcv::PSet>("PixelAnalysis");
  float fthreshold   = pixana_cfg.get<float>("PixelThreshold");
  int fvertex_radius = pixana_cfg.get<int>("PixelRadius");
  int verbosity      = pixana_cfg.get<int>("Verbosity",0);


	return 0;
}