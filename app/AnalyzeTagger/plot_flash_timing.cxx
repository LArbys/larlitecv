#include <iostream>
#include <vector>
#include <sstream>

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/trigger.h"
#include "DataFormat/track.h"
#include "DataFormat/user_info.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "ContainedROI/TaggerFlashMatchAlgoConfig.h"
#include "ContainedROI/TaggerFlashMatchAlgo.h"
#include "ContainedROI/FlashMatchMetricMethods.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "MCTruthTools/extractTruthMethods.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TVector.h"

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
  //enum Stages_t { kThruMu=0, kStopMu, kUntagged, kCROI, kNumStages };
  int kThruMu, kStopMu, kUntagged, kCROI, kNumStages;
  bool use_streclustered = pset.get<bool>("UseReclusteredStopThru");
  std::vector<std::string> stages_track_producers;
  if ( use_streclustered ) {
    kNumStages = 4;
    kThruMu = 0;
    kStopMu = 1;    
    kUntagged = 2;
    kCROI = 3;
    stages_track_producers.resize(kNumStages);
    stages_track_producers[0] = "mergedthrumu3d";
    stages_track_producers[1] = "mergedstopmu3d";
    stages_track_producers[2] = "mergeduntagged3d";
    stages_track_producers[3] = "croi3d";    
  }
  
  // program configuration parameters
  std::vector<int> beam_tick_range(2);
  beam_tick_range[0] = 100;
  beam_tick_range[1] = 400.0;
  const float us_per_tick = 0.015625;
  const float bbox_buffer = 0.0;
  const float match_radius = 20.0;

  // space charge correction class
  larlitecv::SpaceChargeMicroBooNE sce;

  // set up an algo
  larcv::PSet flashalgo_pset = pset.get<larcv::PSet>("TaggerFlashMatchAlgo");
  larlitecv::TaggerFlashMatchAlgoConfig flashalgo_cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( flashalgo_pset );
  larlitecv::TaggerFlashMatchAlgo flashalgo( flashalgo_cfg );
  
  // Flash Match Metrics we want to test

  TFile* rfile = new TFile( pset.get<std::string>("OutputAnaFile").c_str(), "recreate");
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
  int num_true_bbox = 0;
  int num_vtxmatched_tracks = 0;
  float closest_vtx_dist = -1;
  tree->Branch("num_true_bbox", &num_true_bbox, "num_true_bbox/I" );
  tree->Branch("num_vtxmatched_tracks", &num_vtxmatched_tracks, "num_vtxmatched_tracks/I" );
  tree->Branch("closest_vtx_dist", &closest_vtx_dist, "closest_vtx_dist/F" );    
  
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

  // scanned chi2 for minimum
  std::vector<float> smallest_scan_falseflash_v; // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<float> smallest_scan_trueflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  tree->Branch("smallest_scan_falseflashes", &smallest_scan_falseflash_v );
  tree->Branch("smallest_scan_trueflashes",  &smallest_scan_trueflash_v );
  
  // obvious containment
  std::vector<int> containment_trueflash_v;    // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<int> containment_falseflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  tree->Branch("containment_falseflashes", &containment_falseflash_v );
  tree->Branch("containment_trueflashes",  &containment_trueflash_v );

  // cosmic ratio
  std::vector<int> cosmicratio_trueflash_v;    // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<int> cosmicratio_falseflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  tree->Branch("cosmicratio_falseflashes", &cosmicratio_falseflash_v );
  tree->Branch("cosmicratio_trueflashes",  &cosmicratio_trueflash_v );
  

  // containment metric: most negative dwall value, from bounding box
  std::vector<float> roi_dwallx;
  std::vector<float> roi_dwally;
  std::vector<float> roi_dwallz;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;
  truthdata.bindToTree( tree );

  int nentries = dataco[kCROIfile].get_nentries("larcv");

  std::cout << "Start event loop with [enter]" << std::endl;
  std::cin.get();

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
    smallest_gausll_trueflash_v.clear();
    smallest_gausll_falseflash_v.clear();
    smallest_scan_trueflash_v.clear();
    smallest_scan_falseflash_v.clear();    
    containment_trueflash_v.clear();
    containment_falseflash_v.clear();
    cosmicratio_trueflash_v.clear();
    cosmicratio_falseflash_v.clear();
    num_vtxmatched_tracks = 0;
    num_true_bbox = 0;
    closest_vtx_dist = -1;

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

    // get cut data results
    // ---------------------
    larlite::event_user* ev_user_info  = (larlite::event_user*)dataco[kCROIfile].get_larlite_data( larlite::data::kUserInfo, "croicutresults" );    

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
        totpe_hypo += flash.PE(ich)*flashalgo_cfg.fudge_factor;
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
      std::cout << "Loding stage #" << istage << ": " << stages_track_producers[istage] << std::endl;
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
      std::stringstream ss;
      ss << "Number of tracks, " << track_list.size() << ", and number of flash hypotheses, " << flash_hypo_v.size() << ", do not match." << std::endl;
      throw std::runtime_error( ss.str() );
    }
    else if ( flash_hypo_v.size()==0 )
      continue;

    // loop through tracks and flash hypotheses. Get metrics for each.
    for ( size_t itrack=0; itrack<track_list.size(); itrack++ ) {
      const larlite::track& track = *(track_list.at(itrack));
      const larlite::opflash& ophypo = ev_data_ophypo->at(itrack);
      float totpe_data = 0;
      float totpe_hypo = 0.;
      float smallest_chi2 = larlitecv::CalculateFlashMatchChi2( intime_data, ophypo, totpe_data, totpe_hypo, flashalgo_cfg.fudge_factor, false );

      float tot_temp = 0;
      float tot_temp2 = 0;
      float smallest_gausll = larlitecv::CalculateShapeOnlyUnbinnedLL( intime_data, ophypo, tot_temp, tot_temp2, false );

      float smallest_scanchi2 = 0;
      if ( pset.get<bool>("RunScanChi2") )
	smallest_scanchi2 = larlitecv::ScanFlashMatchChi2( 25.0, 25.0, 25.0, 10.0, intime_data, track, flashalgo );

      // -----------------------------------------------
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

      // -----------------------------------------------
      // .. scan tracks for proximinty to neutrino vertex
      bool matchvtx = false;
      int npts = track.NumberTrajectoryPoints();      
      
      for (int ipt=0; ipt<npts; ipt++) {
	float dist = 0.;
	const TVector3& pos = track.LocationAtPoint(ipt);
	for (int v=0; v<3; v++)
	  dist += (pos[v]-apparent_vtx[v])*(pos[v]-apparent_vtx[v]);
	dist = sqrt(dist);
	if ( closest_vtx_dist<0 || closest_vtx_dist>dist )
	  closest_vtx_dist = dist;
	if ( dist<match_radius )
	  matchvtx = true;
      }
      
      if ( matchvtx ) {
	num_vtxmatched_tracks++;
      }
      
      // ------------------------------------
      // Containment cut
      bool contained = true;
      if ( bbox[0][0]<flashalgo_cfg.FVCutX[0] || bbox[0][1]>flashalgo_cfg.FVCutX[1] )
	contained = false;
      if ( bbox[1][0]<flashalgo_cfg.FVCutY[0] || bbox[1][1]>flashalgo_cfg.FVCutY[1] )
	contained = false;
      if ( bbox[2][0]<flashalgo_cfg.FVCutZ[0] || bbox[2][1]>flashalgo_cfg.FVCutZ[1] )
	contained = false;

      // ------------------------------------
      // Cosmic Ratio Cut
      bool cosmicratio = false;
      if ( ev_user_info->size()>0 && (*(ev_user_info->front().get_iarray( "cosmicratio" )))[itrack]==1 )
	cosmicratio = true;

      // ------------------------------------
      // Save Flash Info: True and False
      
      if ( matchvtx ) {
        smallest_chi2_trueflash_v.push_back( smallest_chi2 );
        smallest_gausll_trueflash_v.push_back( smallest_gausll );
	smallest_scan_trueflash_v.push_back( smallest_scanchi2 );
        totpe_hypo_trueflash_v.push_back( totpe_hypo );
        num_true_bbox++;
	if ( contained )
	  containment_trueflash_v.push_back( 1 );
	else
	  containment_trueflash_v.push_back( 0 );

      	if ( cosmicratio )
	  cosmicratio_trueflash_v.push_back( 1 );
	else
	  cosmicratio_trueflash_v.push_back( 0 );
      }
      else {
        smallest_chi2_falseflash_v.push_back( smallest_chi2 );
        smallest_gausll_falseflash_v.push_back( smallest_gausll );
	smallest_scan_falseflash_v.push_back( smallest_scanchi2 );	
        totpe_hypo_falseflash_v.push_back( totpe_hypo );
	if ( contained )
	  containment_falseflash_v.push_back( 1 );
	else
	  containment_falseflash_v.push_back( 0 );	

      	if ( cosmicratio )
	  cosmicratio_falseflash_v.push_back( 1 );
	else
	  cosmicratio_falseflash_v.push_back( 0 );	
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
