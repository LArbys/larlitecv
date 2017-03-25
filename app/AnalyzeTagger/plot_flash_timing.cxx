#include "Base/DataCoordinator.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/trigger.h"
#include "LArUtil/Geometry.h"

#include "TFile.h"
#include "TTree.h"



int main( int nargs, char** argv ) {

  // run pixel analysis. use 

  std::string larlite_input = argv[1];
  

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( larlite_input, "larlite" ); // segment image/original image
  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_source.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larlite") << std::endl;
  
  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );

  TFile* rfile = new TFile("output_flash_timing.root", "recreate");
  TTree* tree = new TTree("flashtimes", "Flash Times");
  int run, subrun, event;
  int flash_time_tick;
  float flash_time_us;  // dwll for end of lepton
  float flash_dtrigger_us;
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");
  tree->Branch("flash_time_tick", &flash_time_tick, "flash_time_tick/I" );
  tree->Branch("flash_time_us", &flash_time_us, "flash_time_us/F" );
  tree->Branch("flash_dtrigger_us", &flash_dtrigger_us, "flash_dtrigger_us/F" );

  int nentries = dataco_source.get_nentries("larlite");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco_source.goto_entry(ientry,"larlite");

    dataco_source.get_id(run,subrun,event);

    if ( ientry%10==0 ) {
      std::cout << "entry " << ientry << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    // get other information, e.g. truth
    larlite::event_opflash* ev_opflash = (larlite::event_opflash*)dataco_source.get_larlite_data(larlite::data::kOpFlash,"simpleFlashBeam");
    larlite::trigger* ev_trigger = (larlite::trigger*)dataco_source.get_larlite_data(larlite::data::kTrigger, "triggersim" );
    for ( auto const& opflash : *ev_opflash ) {
      flash_time_us = opflash.Time();
      flash_dtrigger_us = flash_time_us - ev_trigger->TriggerTime()*0.001;
      flash_time_tick = (int)(flash_time_us/0.015625);
      tree->Fill();
    }
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
