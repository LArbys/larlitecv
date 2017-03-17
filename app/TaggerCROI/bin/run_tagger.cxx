#include <iostream>
#include <string>

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

int main(int nargs, char** argv ) {

  std::cout << "Cosmic Muon Tagger and Contained ROI Selection" << std::endl;

  if ( nargs!=2 ) {
  	std::cout << "usage: ./run_tagger [config file]" << std::endl;
  	return 0;
  }

  std::string cfg_file = argv[1];

  // configuration
  larcv::PSet cfg  = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg.get<larcv::PSet>("TaggerCROI");
  larlitecv::TaggerCROIAlgoConfig tagger_cfg = larlitecv::TaggerCROIAlgoConfig::makeConfigFromFile( cfg_file );


  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  dataco.set_filelist( pset.get<std::string>("LArCVInputFilelist"),   "larcv"   );
  dataco.set_filelist( pset.get<std::string>("LArLiteInputFilelist"), "lalrite" );

  // configure
  dataco.configure( cfg_file, "StorageManager", "IOManager", "TaggerCROI" );
  
  // initialize
  dataco.initialize();

  // Get run entries
  int nentries = dataco.get_nentries("larcv");
  int user_nentries =   pset.get<int>("NumEntries",-1);
  int user_startentry = pset.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }


  larlitecv::TaggerCROIAlgo tagger_algo( tagger_cfg );



  return 0;

}
