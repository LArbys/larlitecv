#include <iostream>
#include <vector>

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/Image2D.h"
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
#include "ContainedROI/TaggerFlashMatchAlgoConfig.h"
#include "ContainedROI/FlashMatchMetricMethods.h"
#include "ContainedROI/FlashROIMatching.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "extractTruthMethods.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

int main( int nargs, char** argv ) {

  std::cout << "STUDY FLASH-MATCH METRICS" << std::endl;

  if ( nargs!=3 ) {
    std::cout << "usage: ./plot_flash_timing [config file #1] [config file #2]" << std::endl;
    return 0;
  }

  // Print Statement #0.5
  std::cout << "Entering the first config file." << std::endl;

  // Load Config File
  // The file that I will want to enter here is: 'flashmetrics.cfg'.
  std::string cfg_file = argv[1];
  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("FlashMatchAnalysis");

  // Print Statement #0.75
  std::cout << "Entering the section of code with the data coordinators." << std::endl;

  // CRB Code.
  // Do the same for the other config file.
  // The file that I want to enter here is: '../ContainedROI/bin/containedroi.cfg'.
  std::string   FlashMgr_cfg_file      = argv[2];
  // Print Statement #0.9
  std::cout << "Entering the statement with the PSetConfiguration." << std::endl;
  larcv::PSet   FlashMgr_cfg_root      = larcv::CreatePSetFromFile( FlashMgr_cfg_file );
  larcv::PSet   FlashMgr_PSet          = FlashMgr_cfg_root.get<larcv::PSet>("ContainedROI").get<larcv::PSet>("TaggerFlashMatchAlgo");
  larlitecv::TaggerFlashMatchAlgoConfig cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( FlashMgr_PSet ); // Give to the instance of class 'TaggerFlashMatchAlgo'.              

  // Declare an instance of class 'TaggerFlashMatchAlgo'.  This is the object that I will use to generate a new flash hypothesis for each of the 'Q_Clusters'.
  larlitecv::TaggerFlashMatchAlgo flash_matcher_algo( cfg );

  // Declare a constant at the top for the number of cm in the x-direction per tick in the wire-plane images.
  // const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec] * [usec/tick] = [cm/tick]

  // Print Statement #1
  std::cout << "Making it past the PSet initialization." << std::endl;

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

  // Print Statement #2
  std::cout << "Making it past the input file list initialization." << std::endl;
  
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
  const float bbox_buffer = 0.0;

  // CRB: Include the configuration parameters for the region that I will check for the lowest chisquare value:
  double LowestChiSquareGridHalfSize     = pset.get< double >("LowestChiSquareGridHalfSize");
  double LowestChiSquareGridPointSep     = pset.get< double >("LowestChiSquareGridPointSep");

  // Print Statement #3
  std::cout << "Making it past the stages track producers initialization." << std::endl;

  // Declare variables for the number of columns and the numbers of rows.
  // size_t row_num = 3456;
  // size_t col_num = 1008;

  // CRB: Declare an 'Image2D' object to reference for generating rows and columns from the wire number and tick.
  // const larcv::Image2D img(row_num, col_num);

  // space charge correction class
  larlitecv::SpaceChargeMicroBooNE sce;

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

  // I will declare a new output file for all of my data.
  TFile* rfile_ana = new TFile("flash_hypo_output_500_events_ana_truth_experiment_100_events.root", "recreate");
  TTree* tree_ev   = new TTree("evtree", "Looking at the Data Event By Event");
  TTree* tree_ana  = new TTree("anatree", "Analysis of the Tagger's Performance");

  // Variables for the data that I will place into my tree.
  std::vector <double> chi_squared_value;
  std::vector <double> x_axis_offset, y_axis_offset, z_axis_offset;
  // double x_vertex_value_sce_and_offset_corrected;
  double x_vertex_value;
  double x_vertex_point_check;
  double x_point_value;
  double y_point_value;
  double z_point_value; // These both have to be declared here so that they can be used in the truth declaration further down.
  double min_chi2_value;
  int    truth_vtx;
  int    event_number;
  int    track_number;
  
  // Declare values for the branches of the 'tree_ev'.
  tree_ev->Branch("event_number", &event_number, "event_number/I");
  tree_ev->Branch("track_number", &track_number, "track_number/I");
  tree_ev->Branch("chi_squared_value", &chi_squared_value);
  tree_ev->Branch("x_axis_offset", &x_axis_offset);
  tree_ev->Branch("y_axis_offset", &y_axis_offset);
  tree_ev->Branch("z_axis_offset", &z_axis_offset);

  // Declare values for the branches of the 'tree_ana'.
  // tree_ana->Branch("x_vertex_value_sce_and_offset_corrected", &x_vertex_value_sce_and_offset_corrected, "x_vertex_value_sce_and_offset_corrected/F");
  tree_ana->Branch("x_vertex_value", &x_vertex_value, "x_vertex_value/D");
  tree_ana->Branch("x_vertex_point_check", &x_vertex_point_check, "x_vertex_point_check/D"); // This is a check comparing the vertex to the first location in the track.
  tree_ana->Branch("min_chi2_value", &min_chi2_value, "min_chi2_value/D");
  tree_ana->Branch("truth_vtx", &truth_vtx, "truth_vtx/I");

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

  // obvious containment
  std::vector<int> containment_trueflash_v;    // chi2 to flash hypothesis not corresponding to true source of trigger
  std::vector<int> containment_falseflash_v;  // chi2 to flash hypothesis corresponding to true source of trigger
  tree->Branch("containment_falseflashes", &containment_falseflash_v );
  tree->Branch("containment_trueflashes",  &containment_trueflash_v );


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

  // Import the larlite

  for (int ientry=0; ientry<nentries; ientry++) {

    // Set the value of 'event_number'.
    event_number = ientry;

    std::cout << "Entering into the event loop for event #" << ientry << "." << std::endl;
    std::cout << "\n" << std::endl;

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
    containment_trueflash_v.clear();
    containment_falseflash_v.clear();

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
    std::cout << "Number of in-time flash hypotheses: " << flash_hypo_v.size() << ". Num zero-pe hypotheses=" << nflash_cut << std::endl;

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

      // Set the value of 'track_number'.
      track_number = itrack;

      // Clear out the vectors that contain the offset information.
      chi_squared_value.clear();
      x_axis_offset.clear();
      y_axis_offset.clear();
      z_axis_offset.clear();
      
      const larlite::track& track = *(track_list.at(itrack));
      const larlite::opflash& ophypo = ev_data_ophypo->at(itrack);
      float totpe_data = 0;
      float totpe_hypo = 0.;
      float smallest_chi2 = larlitecv::CalculateFlashMatchChi2( intime_data, ophypo, totpe_data, totpe_hypo, 20000.0, false );

      float tot_temp = 0;
      float tot_temp2 = 0;
      float smallest_gausll = larlitecv::CalculateShapeOnlyUnbinnedLL( intime_data, ophypo, tot_temp, tot_temp2, false );

      // Implement the code for changing the chi2 test at this point in the code.

      // First, declare a vector for the offset points at the new point in the square around the central point on the track.
      // std::vector < TVector3 > offset_points_vector;
      // std::vector < double > y_offset_vector;
      // std::vector < double > z_offset_vector;

      // Declare a vector for the chi2 values at the different points on the track's trajectory.
      std::vector < float > chi_squared_v;

      // Declare a vector for the number of iterations over the points on the grid.
      int num_iterations = 0;

      larlite::track shifted_track(track);
      
      // The points that are being referenced are with respect to the initial location of these points on the track.
      // Check to see if the change in chi square is constant at all points on the grid.
      for (double x_offset = -LowestChiSquareGridHalfSize; x_offset <= LowestChiSquareGridHalfSize; x_offset += LowestChiSquareGridPointSep) {

	for (double y_offset = -LowestChiSquareGridHalfSize; y_offset <= LowestChiSquareGridHalfSize; y_offset += LowestChiSquareGridPointSep) {

	  for (double z_offset = -LowestChiSquareGridHalfSize; z_offset <= LowestChiSquareGridHalfSize; z_offset += LowestChiSquareGridPointSep) {

	    // Initialize a new track object with a set of trajectory points as given by the 'offset_vector_points'.                                                               
	    // larlite::track shifted_track(track);
	    shifted_track.clear_data();
	    // size_t num_of_trajectory_points = track.NumberTrajectoryPoints();
	    // shifted_track.reserve(num_of_trajectory_points);

	    // Declare an object of type 'TVector3' that will contain the shifted points on the track which will then be appended to 'offset_points_vector'.
	    TVector3 shifted_points_set;

	    // Declare a boolean value for if the track is out of range.
	    bool isShiftedTrackOutOfRange = false;

	    // Fill the vectors above with the offsets of the points of the current track.
	    for (unsigned int point_iterator = 0; point_iterator < track.NumberTrajectoryPoints(); point_iterator++) {

	      // Check to see if any point on the track is out of range.  If it is, then we will skip this particular (x,y,z) shift of the track.
	      // I will put some allowance on the track in the x-direction by allowing it to be placed anywhere from (-50 cm, 300 cm) because the tracks
	      // are not t0-tagged.
	      if ((track.LocationAtPoint(point_iterator)[0] + x_offset) < -50 || (track.LocationAtPoint(point_iterator)[0] + x_offset) > 300 || (track.LocationAtPoint(point_iterator)[1] + y_offset) < -116.5 || (track.LocationAtPoint(point_iterator)[1] + y_offset) > 116.5 || (track.LocationAtPoint(point_iterator)[2] + z_offset) < 0.0 || (track.LocationAtPoint(point_iterator)[2] + z_offset) > 1036.8) {
		
		// Reset the value of 'isShiftedTrackOutOfRange' to 'true'.
		isShiftedTrackOutOfRange = true;
	       
		break;

	      }

	      // Offset each of the points on the track, append them to 'shifted_points_set', and append that object to 'offset_points_vector'.
	      shifted_points_set[0] = track.LocationAtPoint(point_iterator)[0] + x_offset;
	      shifted_points_set[1] = track.LocationAtPoint(point_iterator)[1] + y_offset;
	      shifted_points_set[2] = track.LocationAtPoint(point_iterator)[2] + z_offset;

	      // Append these points onto the 'shifted_track'.
	      shifted_track.add_vertex(shifted_points_set);
	      shifted_track.add_direction(track.DirectionAtPoint(point_iterator));	      

	      // Append this object onto the end of the 'offset_points_vector' object (this will occur in order so it will not be a problem).
	      // offset_points_vector.push_back(shifted_points_set);

	    }

	    // If 'isShiftedTrackOutOfRange' is now 'true', then you can continue onto the next set of shifted points that we will consider.
	    if (isShiftedTrackOutOfRange)
	      continue;

	    // Turn this track into an object of type 'flashana::QCluster_t' by using the function 'TaggerFlashMatchAlgo::GenerateQCluster'.
	    flashana::QCluster_t shifted_track_qcluster = flash_matcher_algo.GenerateQCluster(shifted_track);

	    // Generate an unfitted hypothesis of type 'flashana::Flash_t' for this 'qcluster' by using the 'TaggerFlashMatchAlgo::GenerateUnfittedFlashHypothesis' function.
	    flashana::Flash_t shifted_track_flash = flash_matcher_algo.GenerateUnfittedFlashHypothesis(shifted_track_qcluster);

	    // Generate this 'flashana::Flash_t' into an object of type 'larlite::opflash' with the function 'TaggerFlashMatchAlgo::MakeOpFlashFromFlash'
	    larlite::opflash shifted_track_ophypo = flash_matcher_algo.MakeOpFlashFromFlash(shifted_track_flash);
	   
	    // Generate a value for the minimum chi squared value with the variales 'totpe_data_shifted' and 'totpe_hypo_shifted' that are used.
	    float totpe_data_shifted = 0.;
	    float totpe_hypo_shifted = 0.;

	    float smallest_chi2_shifted = larlitecv::CalculateFlashMatchChi2( intime_data, shifted_track_ophypo, totpe_data_shifted, totpe_hypo_shifted, 20000.0, false );
	    num_iterations++;

	    // Set the value for the chi-squared calculation.
	    chi_squared_value.push_back(smallest_chi2_shifted);
	    x_axis_offset.push_back(x_offset);
	    y_axis_offset.push_back(y_offset);
	    z_axis_offset.push_back(z_offset);
	    
	    // std::cout << "The new chi squared value = " << smallest_chi2_shifted << "." << std::endl;
	    // std::cout << "The number of iterations this is = " << num_iterations << "." << std::endl;
	    // std::cout << "\n" << std::endl;

	    // Append this value onto 'chi_squared_v'.
	    chi_squared_v.push_back(smallest_chi2_shifted);

	  }

	}

      }

      // Print out 'chi_squared_v'.
      // std::cout << "The length of the vector of chi squared values is: " << chi_squared_v.size() << "." << std::endl;
      // std::cout << "\n" << std::endl;     

      // I will run a loop to make sure the output that we have is minimized. I declared this variable further up in my code as well to avoid declaration errors.
      min_chi2_value  = 10000.;
      // TVector3 min_chi2_coords = shifted_track.LocationAtPoint(0);

      // Loop through the chi2 values to see what the minimum value is.
      for (int chi2_iter = 0; chi2_iter < chi_squared_v.size(); chi2_iter++) {

	// std::cout << "chi2 iter: " << chi2_iter << "." << std::endl;
	// std::cout << "chi2 value: " << chi_squared_v.at(chi2_iter) << "." << std::endl;
	// std::cout << "\n" << std::endl;

	// Check to see if the value of 'chi_squared_v' at this index is less than 'min_chi2_value'.  Set it equal to 'min_chi2_value' if it is.
	if (chi_squared_v.at(chi2_iter) < min_chi2_value) {

	  // std::cout << "The new minimum chi2 value = " << chi_squared_v.at(chi2_iter) << "." << std::endl;

	  min_chi2_value  = chi_squared_v.at(chi2_iter);
	  // min_chi2_coords = shifted_track.LocationAtPoint(chi2_iter);

	}

      }

      // std::cout << "The minimum chi2 value in the vicinity of this track is = " << min_chi2_value << "." << std::endl;
      // std::cout << "\n" << std::endl;

      // std::cout << "The chi2 value for this track at its original location = " << smallest_chi2 << "." << std::endl;
      // std::cout << "\n" << std::endl;

	    // Make a copy of the current value of the track using the copy track constructor.
            // track(const track& original) :   
            //  data_base(original),
	    //  fXYZ(original.fXYZ),
	    //  fDir(original.fDir),
	    //  fCov(original.fCov),
	    //  fdQdx(original.fdQdx),
	    //  fFitMomentum(original.fFitMomentum),
	    //  fID(original.fID)
      
      // need to determine if this track is the intime track or not
      // we do this by checking if track gets close to true vertex. slow AF.
      std::vector< std::vector<float> > bbox = larlitecv::TaggerFlashMatchAlgo::GetAABoundingBox( track );
      bool inbbox = true;
      std::vector<double> sce_offsets = sce.GetPosOffsets( truthdata.pos[0], truthdata.pos[1], truthdata.pos[2] );
      std::vector<double> apparent_vtx(3,0);
      apparent_vtx[0] = truthdata.pos[0] - sce_offsets[0] - 0.17;
      apparent_vtx[1] = truthdata.pos[1] + sce_offsets[1];
      apparent_vtx[2] = truthdata.pos[2] + sce_offsets[2];

      // I can set the values for the truth vertex and the 'corrected' vertex here.
      // x_vertex_value                          = truthdata.pos[0];
      // x_vertex_value_sce_and_offset_corrected = apparent_vtx[0];

      x_vertex_value       = track.Vertex()[0];
      x_vertex_point_check = track.LocationAtPoint(0)[0];
      // y_vertex_value = track.Vertex()[1];
      // z_vertex_value = track.Vertex()[2];
      

      // Print these two values.
      // std::cout << "x vertex value = " << x_vertex_value << "." << std::endl;
      // std::cout << "x corrected vertex value = " << x_vertex_value_sce_and_offset_corrected << "." << std::endl;

      // Set 'truth_vtx' to 0.
      truth_vtx = 0;

      for (size_t track_point_iter = 0; track_point_iter < track.NumberTrajectoryPoints(); track_point_iter++) {

	// Declare the points that will be compared to the 'apparent_vtx' values.
	x_point_value = track.LocationAtPoint(track_point_iter)[0];
	y_point_value = track.LocationAtPoint(track_point_iter)[1];
	z_point_value = track.LocationAtPoint(track_point_iter)[2];

	// Check if this is a true vertex by seeing if it is 5 cm from the 'apparent_vtx' values.
	if (fabs(apparent_vtx[0] - x_point_value) < 5.0 && fabs(apparent_vtx[1] - y_point_value) < 5.0 && fabs(apparent_vtx[2] - z_point_value) < 5.0)
	  truth_vtx = 1;

      }
      
      for ( int v=0; v<3; v++ ) {
        if ( apparent_vtx[v]<bbox[v][0]-bbox_buffer || apparent_vtx[v]>bbox[v][1]+bbox_buffer ) {
          inbbox = false;
          break;
        }
      }

      bool contained = false;
      if ( bbox[0][0]>=-10 && bbox[0][0]<=275 && bbox[0][1]>=-10 && bbox[0][1]<=275 )
	contained = true;

      if ( inbbox ) {

	// If this loop is entered, then it is a real vertex. Set 'truth_vtx' to 1.
	// truth_vtx = 1;
	
        smallest_chi2_trueflash_v.push_back( smallest_chi2 );
        smallest_gausll_trueflash_v.push_back( smallest_gausll );
        totpe_hypo_trueflash_v.push_back( totpe_hypo );
        num_true_bbox++;
	if ( contained )
	  containment_trueflash_v.push_back( 1 );
	else
	  containment_trueflash_v.push_back( 0 );
      }
      else {
        smallest_chi2_falseflash_v.push_back( smallest_chi2 );
        smallest_gausll_falseflash_v.push_back( smallest_gausll );
        totpe_hypo_falseflash_v.push_back( totpe_hypo );
	if ( contained )
	  containment_falseflash_v.push_back( 1 );
	else
	  containment_falseflash_v.push_back( 0 );	
      }

      // Print the x-coordinate value.
      std::cout << "x-coordinate value = " << x_vertex_value << "." << std::endl;
      std::cout << "minimum chi2 value = " << min_chi2_value << "." << std::endl;
      std::cout << "\n" << std::endl;

      // Fill the 'tree_ana' and 'tree_ev' here for each track in the set of events.
      tree_ana->Fill();
      tree_ev->Fill();
      
    }
    std::cout << "Number of Vertex Enclosing BBox: " << num_true_bbox << std::endl;

    tree->Fill();

    // I will only use 10 entries for now.
    if ( ientry>=99)
      break;
  }//end of entry loop

  rfile->Write();
  rfile_ana->Write();
  return 0;

}//end of main
