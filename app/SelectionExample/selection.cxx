#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <ctime>

// config/storage: from LArCV
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h" // example of data product
#include "DataFormat/mctruth.h" // example of data product
#include "DataFormat/mcnu.h" // example of data product
#include "DataFormat/mcpart.h" // example of data product
#include "DataFormat/mctrack.h" // example of data product

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"



int main( int nargs, char** argv ) {
  
  std::cout << "[Example LArCV and LArLite Data Access with Event Selection]" << std::endl;
  
  // configuration
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet select_config = cfg.get<larcv::PSet>("SelectionConfigurationFile");
  std::string larcv_image_producer = select_config.get<std::string>("InputLArCVImages");
  std::string larlite_mctruth_producer = select_config.get<std::string>("InputMCTruthProducer");
  std::string larcv_flist   = select_config.get<std::string>("LArCVFilelist");
  std::string larlite_flist = select_config.get<std::string>("LArLiteFilelist");
  std::vector<int> modes = select_config.get<std::vector<int>>("SelectedModes");
  std::vector<int> currents = select_config.get<std::vector<int>>("SelectedCurrents");
  std::vector<float> Enu_bounds_GeV = select_config.get<std::vector<float>>("EnuBoundsGeV");
  int start_entry = select_config.get<int>("StartEntry", 0);
  int max_entries = select_config.get<int>("MaxEntries",-1);

  if ( Enu_bounds_GeV.size()!=2 ) {
    throw std::runtime_error("EnuBounds_GeV must have two values.");
  }
  
  // Configure Data coordinator
  larlitecv::DataCoordinator dataco;

  // you can get example data files from: /uboone/data/users/tmw/tutorials/larlitecv_example/

  // ----------------------------------------------------------------------------------------------------------------
  // // individual files, for debug
  // // larlite
  // //dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/larlite/larlite_wire_0000.root", "larlite" );
  // dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/larlite/larlite_chstatus_0000.root", "larlite" );
  // dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/larlite/larlite_opdigit_0000.root", "larlite" );
  // dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/larlite/larlite_opreco_0000.root", "larlite" );
  // dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/larlite/larlite_mcinfo_0000.root", "larlite" );

  // // larcv
  // dataco.add_inputfile( "/a/data/amsterdam/tmw23/sw/../condor/thrumu/mcc7_bnb_cosmic_v00_p00b/supera/supera_0000.root", "larcv" );
  // ----------------------------------------------------------------------------------------------------------------
  // use filelist
  dataco.set_filelist( larcv_flist,   "larcv" );
  dataco.set_filelist( larlite_flist, "larlite" );

  // configure
  dataco.configure( "config.cfg", "StorageManager", "IOManager", "SelectionConfigurationFile" );
  
  // initialize
  dataco.initialize();


  // Start Event Loop
  int nentries = dataco.get_nentries("larcv");
  int end_entry = nentries;
  if ( max_entries>=0 ) {
    end_entry = start_entry + max_entries;
    end_entry = ( end_entry > nentries ) ? nentries : end_entry;
  }

  for (int ientry=start_entry; ientry<end_entry; ientry++) {
    std::cout << "[Entry " << ientry << "]" << std::endl;

    dataco.goto_entry(ientry,"larcv");

    // we need the MCTruth information
    larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)dataco.get_larlite_data( larlite::data::kMCTruth, larlite_mctruth_producer );
    const larlite::mcnu& neutrino = ev_mctruth->at(0).GetNeutrino();

    // get the larcv data
    larcv::EventImage2D* ev_img2d = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "tpc" );
    const std::vector<larcv::Image2D>& img_v = ev_img2d->Image2DArray();

    // is it in the modes?
    bool modefound = false;
    for ( auto& mode : modes ) {
      if ( neutrino.InteractionType()==mode ) {
	modefound = true;
	break;
      }
    }

    bool currentfound = false;
    for (auto& current : currents ) {
      if ( current==neutrino.CCNC() ) {
	currentfound = true;
	break;
      }
    }

    bool withinenergy = false;
    if ( neutrino.Nu().Momentum(0).E()>Enu_bounds_GeV[0] && neutrino.Nu().Momentum(0).E()<Enu_bounds_GeV[1] ) {
      withinenergy = true;
    }

    const TLorentzVector& nu_pos = neutrino.Nu().Position();
    bool withinTPC = false;
    if ( nu_pos.X()>-1.0 && nu_pos.X()<270.0 && nu_pos.Y()>-118.0 && nu_pos.Y()<118.0 && nu_pos.Z()>0 && nu_pos.Z()<1037.0 )
      withinTPC = true;

    bool passes = modefound & currentfound & withinenergy & withinTPC;

    std::cout << "larlite RSE=(" << ev_mctruth->run() << "," << ev_mctruth->subrun() << "," << ev_mctruth->event_id() << ")" << std::endl;
    std::cout << "larcv RSE=(" << ev_img2d->run() << "," << ev_img2d->subrun() << "," << ev_img2d->event() << ")" << std::endl;

    std::cout << "neutrino passes: " << passes << std::endl;
    std::cout << " energy=" << neutrino.Nu().Momentum(0).E() << " pos=(" << nu_pos.X() << "," << nu_pos.Y() << "," << nu_pos.Z() << ")" << std::endl;
    std::cout << " mode=" << neutrino.InteractionType() << " current=" << neutrino.CCNC() << std::endl;


    // go to tree
    if ( passes )
      dataco.save_entry();
  }

  std::cout << "finalize." << std::endl;

  // finalize output files
  dataco.finalize();

  return 0;
}
