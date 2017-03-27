#include "Base/DataCoordinator.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/trigger.h"
#include "DataFormat/track.h"
#include "LArUtil/Geometry.h"

#include "TFile.h"
#include "TTree.h"

int main( int nargs, char** argv ) {

  std::cout << "STUDY FLASH-MATCH METRICS" << std::endl;

  if ( nargs!=2 ) {
    std::cout << "usage: ./plot_flash_timing [config file]" << std::endl;
    return 0;
  }

  // Load Config File
  std::string cfg_file = argv[1];
  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("FlashMatchAnalysis");

  // Setup the Data Coordinators
  enum { kSource=0, kCROIfile, kNumSourceTypes } SourceTypes_t;
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "FlashMatchAnalysis" );
    dataco[isrc].initialize();
  }
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }

  // Stages
  enum { kThruMu=0, kStopMu, kUntagged, kCROI, kNumStages } Stages_t;
  std::string stages_track_producers[kNumStages] = { "thrumu3d", "stopmu3d", "untagged3d", "croi3d" };

  // program configuration parameters
  std::vector<int> beam_tick_range(2);
  beam_tick_range[0] = 100;
  beam_tick_range[1] = 400.0;
  const float us_per_tick = 0.015625;

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

  int nentries = dataco[kCROIfile].get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(run,subrun,event);
    for ( int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
      if ( isrc!=kCROIfile ) {
        dataco[isrc].goto_event(run,subrun,event, "larcv");
      }
    }

    std::cout << "======================================================================" << std::endl;
    std::cout << "entry " << ientry << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    // we get wnat to get:
    //  (1) truth info
    //  (2) "data" in time opflashes
    //  (3) flash hypotheses from tracks
    //  (4) calculate metrics for flash matching
    //  (5) split flash hypotheses coming from "correct" and "incorrect" hypotheses

    // get flash data
    // --------------
    larlite::event_opflash* ev_data_opflash = (larlite::event_opflash*)dataco[kSource].get_larlite_data( larlite::data::kOpFlash, "simpleFlashBeam" );
    larlite::event_opflash* ev_data_ophypo  = (larlite::event_opflash*)dataco[kCROIfile].get_larlite_data( larlite::data::kOpFlash, "ophypo" );

    // Get in time beam flashes
    std::vector<larlite::opflash> intime_data;
    for ( auto const& flash : *ev_data_opflash ) {
      int tick = flash.Time()/us_per_tick;
      if ( tick>=beam_tick_range[0] && tick<=beam_tick_range[1] )
        intime_data.push_back( flash );
    }
    std::cout << "Number of in-time (data) flashes: " << intime_data.size() << std::endl;
    std::cout << "Number of in-time flash hypotheses: " << ev_data_ophypo->size() << std::endl;

    // get track data
    // --------------
    int ntracks = 0;
    larlite::event_track* ev_tagger_tracks[kNumStages];
    std::vector<const larlite::track*> track_list;
    for ( int istage=0; istage<=kUntagged; istage++ ) {
      ev_tagger_tracks[istage] = (larlite::event_track*)dataco[kCROIfile].get_larlite_data( larlite::data::kTrack, stages_track_producers[istage] );
      for ( auto const& track : *(ev_tagger_tracks[istage]) ) {
        if ( track.NumberTrajectoryPoints()==0 )
          continue;
        track_list.push_back( &track );
      }
    }
    std::cout << "Number of tracks: " << track_list.size() << std::endl;

    if ( track_list.size()!=ev_data_ophypo->size() ) {
      throw std::runtime_error("Number of tracks and number of flash hypotheses do not match.");
    }


    if ( ientry>=3)
      break;
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
