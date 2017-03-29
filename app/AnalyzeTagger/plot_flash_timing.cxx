#include <iostream>
#include <vector>

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

// larlitecv
#include "Base/DataCoordinator.h"
#include "ContainedROI/TaggerFlashMatchAlgo.h"
#include "ContainedROI/FlashMatchMetricMethods.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "extractTruthMethods.h"

// ROOT
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
  const float bbox_buffer = 10.0;

  // space charge correction class
  larlitecv::SpaceChargeMicroBooNE sce;

  // Flash Match Metrics we want to test

  TFile* rfile = new TFile("output_flash_timing.root", "recreate");
  TTree* tree = new TTree("flashtimes", "Flash Times");
  int run, subrun, event;

  // details about data flash
  int flash_time_tick;
  float flash_time_us;  // dwll for end of lepton
  float flash_dtrigger_us;
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");
  tree->Branch("flash_time_tick", &flash_time_tick, "flash_time_tick/I" );
  tree->Branch("flash_time_us", &flash_time_us, "flash_time_us/F" );
  tree->Branch("flash_dtrigger_us", &flash_dtrigger_us, "flash_dtrigger_us/F" );


  // flash metrics output

  // binned poisson likelihood ratio
  std::vector<float> smallest_chi2_falseflash_v; // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<float> smallest_chi2_trueflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  std::vector<float> totpe_data_v;               // total pe for in-time flashes in the data
  std::vector<float> totpe_hypo_trueflash_v;     // total pe for correct flash hypotheses
  std::vector<float> totpe_hypo_falseflash_v;    // total pe for flase flash hypotheses
  tree->Branch("smallest_chi2_falseflashes", &smallest_chi2_falseflash_v );
  tree->Branch("smallest_chi2_trueflashes",  &smallest_chi2_trueflash_v );
  tree->Branch("totpe_data",                 &totpe_data_v );
  tree->Branch("totpe_hypo_trueflashes",     &totpe_hypo_trueflash_v );
  tree->Branch("totpe_hypo_falseflashes",    &totpe_hypo_falseflash_v );

  // unbinned negative log-likelihood using multivariate guassian distribution
  std::vector<float> smallest_gausll_falseflash_v; // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<float> smallest_gausll_trueflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  tree->Branch("smallest_gausll_falseflashes", &smallest_gausll_falseflash_v );
  tree->Branch("smallest_gausll_trueflashes",  &smallest_gausll_trueflash_v );


  // containment metric: most negative dwall value, from bounding box
  std::vector<float> roi_dwallx;
  std::vector<float> roi_dwally;
  std::vector<float> roi_dwallz;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;
  truthdata.bindToTree( tree );

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

    // clear variables
    truthdata.clear();
    smallest_chi2_trueflash_v.clear();
    smallest_chi2_falseflash_v.clear();
    totpe_data_v.clear();
    totpe_hypo_trueflash_v.clear();
    totpe_hypo_falseflash_v.clear();

    // we get wnat to get:
    //  (1) truth info
    //  (2) "data" in time opflashes
    //  (3) flash hypotheses from tracks
    //  (4) calculate metrics for flash matching
    //  (5) split flash hypotheses coming from "correct" and "incorrect" hypotheses

    // get truth data
    // --------------
    larlite::event_mctruth* ev_mctruth   = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
    larlite::event_mctrack* ev_mctrack   = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");
    larlite::trigger* ev_trigger         = (larlite::trigger*)       dataco[kSource].get_larlite_data(larlite::data::kTrigger,"triggersim");
    larlitecv::extractTruth( *ev_mctruth, *ev_mctrack, truthdata );

    // get flash data
    // --------------
    larlite::event_opflash* ev_data_opflash = (larlite::event_opflash*)dataco[kSource].get_larlite_data( larlite::data::kOpFlash, "simpleFlashBeam" );
    larlite::event_opflash* ev_data_ophypo  = (larlite::event_opflash*)dataco[kCROIfile].get_larlite_data( larlite::data::kOpFlash, "ophypo" );

    // select in time beam flashes
    std::vector<larlite::opflash> intime_data;
    for ( auto const& flash : *ev_data_opflash ) {
      int tick = flash.Time()/us_per_tick;
      if ( tick>=beam_tick_range[0] && tick<=beam_tick_range[1] ) {
        intime_data.push_back( flash );
        float totpe_data = 0.;
        for (int ich=0; ich<32; ich++)
          totpe_data += flash.PE(ich);
        totpe_data_v.push_back( totpe_data );
      }
    }
    std::cout << "Number of in-time (data) flashes: " << intime_data.size() << std::endl;

    // select non-zero flash hypotheses
    std::vector<larlite::opflash> flash_hypo_v;
    int nflash_cut = 0;
    for ( auto const& flash : *ev_data_ophypo ) {
      float totpe_hypo = 0.;
      for (int ich=0; ich<32; ich++) {
        totpe_hypo += flash.PE(ich)*20000.0;
      }
      flash_hypo_v.push_back( flash );
    }
    std::cout << "Number of in-time flash hypotheses: " <<flash_hypo_v.size() << ". Num zero-pe hypotheses=" << nflash_cut << std::endl;

    // get track data
    // --------------
    int ntracks = 0;
    int ntracks_cut = 0;
    larlite::event_track* ev_tagger_tracks[kNumStages];
    std::vector<const larlite::track*> track_list;
    for ( int istage=0; istage<=kUntagged; istage++ ) {
      ev_tagger_tracks[istage] = (larlite::event_track*)dataco[kCROIfile].get_larlite_data( larlite::data::kTrack, stages_track_producers[istage] );
      for ( auto const& track : *(ev_tagger_tracks[istage]) ) {
        if ( istage==kUntagged && track.NumberTrajectoryPoints()==0 ) {
        //if ( track.NumberTrajectoryPoints()==0 ) {
          ntracks_cut++;
          continue;
        }
        track_list.push_back( &track );
      }
    }
    std::cout << "Number of tracks: " << track_list.size() << ". number removed=" << ntracks_cut << "." << std::endl;

    if ( flash_hypo_v.size()>0 && track_list.size()!=flash_hypo_v.size() ) {
      //continue; // freak out and skip event
      throw std::runtime_error("Number of tracks and number of flash hypotheses do not match.");
    }
    else if ( flash_hypo_v.size()==0 )
      continue;

    // loop through tracks and flash hypotheses. Get metrics for each.
    int num_true_bbox = 0;
    for ( size_t itrack=0; itrack<track_list.size(); itrack++ ) {
      const larlite::track& track = *(track_list.at(itrack));
      const larlite::opflash& ophypo = ev_data_ophypo->at(itrack);
      float totpe_data = 0;
      float totpe_hypo = 0.;
      float smallest_chi2 = larlitecv::CalculateFlashMatchChi2( intime_data, ophypo, totpe_data, totpe_hypo, 20000.0, false );

      float tot_temp = 0;
      float tot_temp2 = 0;
      float smallest_gausll = larlitecv::CalculateShapeOnlyUnbinnedLL( intime_data, ophypo, tot_temp, tot_temp2, false );

      // need to determine if this track is the intime track or not
      // we do this by checking if track gets close to true vertex. slow AF.
      std::vector< std::vector<float> > bbox = larlitecv::TaggerFlashMatchAlgo::GetAABoundingBox( track );
      bool inbbox = true;
      std::vector<double> sce_offsets = sce.GetPosOffsets( truthdata.pos[0], truthdata.pos[1], truthdata.pos[2] );
      std::vector<double> apparent_vtx(3,0);
      apparent_vtx[0] = truthdata.pos[0] - sce_offsets[0] - 0.17;
      apparent_vtx[1] = truthdata.pos[1] + sce_offsets[1];
      apparent_vtx[2] = truthdata.pos[2] + sce_offsets[2];
      for ( int v=0; v<3; v++ ) {
        if ( apparent_vtx[v]<bbox[v][0]-bbox_buffer || apparent_vtx[v]>bbox[v][1]+bbox_buffer ) {
          inbbox = false;
          break;
        }
      }

      if ( inbbox ) {
        smallest_chi2_trueflash_v.push_back( smallest_chi2 );
        smallest_gausll_trueflash_v.push_back( smallest_gausll );
        totpe_hypo_trueflash_v.push_back( totpe_hypo );
        num_true_bbox++;
      }
      else {
        smallest_chi2_falseflash_v.push_back( smallest_chi2 );
        smallest_gausll_falseflash_v.push_back( smallest_gausll );
        totpe_hypo_falseflash_v.push_back( totpe_hypo );
      }
    }
    std::cout << "Number of Vertex Enclosing BBox: " << num_true_bbox << std::endl;

    tree->Fill();

    //if ( ientry>=3)
    //  break;
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
