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

int main( int nargs, char** argv ) {

  std::cout << "STUDY FLASH-MATCH METRICS" << std::endl;

  if ( nargs!=2 ) {
    std::cout << "usage: ./plot_flash_timing [config file]" << std::endl;
    return 0;
  }

  // Load Config File
  // The file that I will want to enter here is: 'flashmetrics.cfg'.
  std::string cfg_file = argv[1];
  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("FlashMatchAnalysis");

  // CRB Code.
  // Do the same for the other config file.
  // The file that I want to enter here is: '../ContainedROI/bin/containedroi.cfg'.
  std::string   FlashMgr_cfg_file      = argv[2];
  larcv::PSet FlashMgr_cfg_root        = larcv::CreatePSetFromFile( FlashMgr_cfg_file );
  const larcv::PSet& FlashMgr_pset     = FlashMgr_cfg_root.get<larcv::PSet>("ContainedROI");

  FlashMgr_cfg_root.get<larcv::PSet>("ContainedROI").get<larcv::PSet>(“TaggerFlashMatchAlgo”);
  TaggerFlashMatchAlgo

  // Create an object of class 'TaggerMatchAlgoConfig'.
  larlitecv::TaggerFlashMatchAlgoConfig first_Config_obj;
  
  // Use this object to produce another object of the same class that will initialize the flash manager.
  larlitecv::TaggerFlashMatchAlgoConfig second_Config_obj = first_Config_obj.MakeTaggerFlashMatchAlgoConfigFromPSet( FlashMgr_pset );

  // With this information, I can feed this second object into an object of class 'TaggerFlashMatchAlgo', which will then configure the 'm_flash_matcher' for me.
  // This is the function that I can also use to run the chi hypotheses with as well.
  larlitecv::TaggerFlashMatchAlgo chi_squared_hypothesis_object ( second_Config_obj );

  // Declare a constant at the top for the number of cm in the x-direction per tick in the wire-plane images.
  const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec] * [usec/tick] = [cm/tick]

  larlitecv::TaggerFlashMatchAlgoConfig cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( set ); // Give to the 

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
  const float bbox_buffer = 0.0;

  // CRB: Include the configuration parameters for the region that I will check for the lowest chisquare value:
  double LowestChiSquareGridHalfSize     = pset.get< double >("LowestChiSquareGridHalfSize");
  double LowestChiSquareGridPointSep     = pset.get< double >("LowestChiSquareGridPointSep");

  // Declare variables for the number of columns and the numbers of rows.
  size_t row_num = 3456;
  size_t col_num = 1008;

  // CRB: Declare an 'Image2D' object to reference for generating rows and columns from the wire number and tick.
  const larcv::Image2D& img(row_num, col_num);

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
      const larlite::track& track = *(track_list.at(itrack));
      const larlite::opflash& ophypo = ev_data_ophypo->at(itrack);
      float totpe_data = 0;
      float totpe_hypo = 0.;
      float smallest_chi2 = larlitecv::CalculateFlashMatchChi2( intime_data, ophypo, totpe_data, totpe_hypo, 20000.0, false );

      float tot_temp = 0;
      float tot_temp2 = 0;
      float smallest_gausll = larlitecv::CalculateShapeOnlyUnbinnedLL( intime_data, ophypo, tot_temp, tot_temp2, false );

      // Implement the code for changing the chi2 test at this point in the code.

      // First, declare vectors for the amount of the offset for the track at each new value of the chi2 in each of the coordinates.
      std::vector < double > x_offset_vector;
      std::vector < double > y_offset_vector;
      std::vector < double > z_offset_vector;
      
      // The points that are being referenced are with respect to the initial location of these points on the track.
      // Check to see if the change in chi square is constant at all points on the grid.
      for (double x_coord = -LowestChiSquareGridHalfSize; x_coord <= LowestChiSquareGridHalfSize; x_coord += LowestChiSquareGridPointSep) {

	for (double y_coord = -LowestChiSquareGridHalfSize; y_coord <= LowestChiSquareGridHalfSize; y_coord += LowestChiSquareGridPointSep) {

	  for (double z_coord = -LowestChiSquareGridHalfSize; z_coord <= LowestChiSquareGridHalfSize; z_coord += LowestChiSquareGridPointSep) {

	    // Declare a boolean for if any part of the track is out of range.
	    bool part_of_track_out_of_range = false;

	    // Declare an empty vector of the new points in each dimension.  I will loop through the points and offset each of them by the amount given above.
	    std::vector < double > offset_x_coordinates;
	    std::vector < double > offset_y_coordinates;
	    std::vector < double > offset_z_coordinates;

	    // Loop through the points on the track and offset each of them by the amount given by the '*dim*_coord', and then append those values to the vectors above.
	    // Coordinates corresponding to the same points are appended at the same position in the vectors, which maintains organization of the output.
	    for (unsigned int offset_iterator = 0; offset_iterator < track.size(); offset_iterator++) {

	      // Check the boundary conditions for if any point on the track is out of range of the TPC.
	      // I allow some offset on the track x endpoints, because the tracks are not t0-corrected and therefore could be reconstructed outside of the TPC.
	      if ((track.LocationAtPoint(offset_iterator)[0] + x_coord) < -50.0 || (track.LocationAtPoint(offset_iterator)[0] + x_coord) > 350 || (track.LocationAtPoint(offset_iterator)[1] + y_coord) < -116.5 || (track.LocationAtPoint(offset_iterator)[1] + y_coord) > 116.5 || (track.LocationAtPoint(offset_iterator)[2] + z_coord) < 0.0 || (track.LocationAtPoint(offset_iterator)[2] + z_coord) > 1036.8) {

		// Set 'part_of_track_out_of_range' to 'true'.
		part_of_track_out_of_range = true;

		// Break the loop.
		break;

	      }

	      offset_x_coordinates.push_back(track.LocationAtPoint(offset_iterator)[0] + x_coord);
	      offset_y_coordinates.push_back(track.LocationAtPoint(offset_iterator)[1] + y_coord);
	      offset_z_coordinates.push_back(track.LocationAtPoint(offset_iterator)[2] + z_coord);

	    }

	    // Continue onto the next set of offsets if part of the track is out of range.  
	    if (part_of_track_out_of_range = true)

	      continue;

	  }
	
	    // Loop through each of the planes, and for each plane loop through each of the points to (1) convert them into a tick and a wire number, (2) convert them into a row and column number, and (3) declare a Pixel2D object with these items that can be appended to the 'Pixel2DCluster' lists above.

	    for (size_t plane = 0; plane < 3; plane++) {

	      // Convert each of these coordinates to wire and tick values, which can then be converted into row and column pixel values.                     
	      // This will all happen in the same loop.  Declare a vector of lists of pixel values that you'll need (there will be a separate list for each of the planes).   
	      larcv::Pixel2DCluster pixel_clusters;

	      for (size_t point_iterator = 0; point_iterator < offset_x_coordinates.size(); point_iterator++) {

		// Declare a 'Double_t' object for the shifted coordinates at this point.
		Double_t xyz[3] = {0};

		// Initialize the values of 'xyz' to the values of the coordinate vectors at this index.
		xyz[0] = offset_x_coordinates[point_iterator];
		xyz[1] = offset_y_coordinates[point_iterator];
		xyz[2] = offset_z_coordinates[point_iterator];

		// Use these coordinates and the 'LArUtil' geometry considerations to find the tick (using the x-coordinate) and the wire number (using the yz-coordinates) for this coordinate.
		float fwire = larutil::Geometry::GetME()->WireCoordinate( xyz , plane );
		fwire       = ( fwire<0 ) ? 0 : fwire; // Ensure that the wire number is not negative.
		float tick  = tick  = xyz[0]/cm_per_tick+3200.0;

		// Round off both of these values.
		fwire       = round(fwire);
		tick        = round(tick);
		
		// Produce row and column values using the image functionality from LArCV.
		int row = img.meta().row(tick);
		int col = img.meta().col(fwire);

		// Convert the row and column into a Pixel2D object.
		// Note that this takes input in the order 'x_, y_' which corresponds to an ordering of the Image2D coordinates of (col, row), opposite to what we have been doing before.
		larcv::Pixel2D pixel(col, row);

		// Append this Pixel2D object onto the 'pixel_clusters' array at the entry at index 'plane'.
		pixel_clusters.push_back(pixel);

	      }

	      // Within the loop over the planes, I can now declare a 'constant' vector for the pixel clusters to give the first function the correct input.
	      // const larcv::Pixel2DCluster& constant_pixel_cluster = pixel_clusters;

	      // Declare an 'extended_cluster' that will be fed as input into the next function.
	      // It is a greater cluster of points around the original cluster that is formed in the process of generating a flash.
	      larcv::Pixel2DCluster& expanded_cluster;

	      // Generate a 'QCluster_t' object with the 'MakeTPCFlashObject' function.
	      flashana::QCluster_t current_qcluster = FlashROIMatching::MakeTPCFlashObject(constant_pixel_cluster, img, expanded_cluster);

	      // With this, I can calculate an unfitted hypothesis for the chi squared using the 'current_qcluster' variable.
	      // Place the qcluster output into the same format as it is in the 'GenerateUnfittedFlashHypothesis' function.
	      const flashana::QCluster_t& constant_qcluster = current_qcluster;

	      // Feed this into the function 'GenerateUnfittedFlashHypothesis'.
	      flashana::Flash_t unfittedhypothesis = TaggerFlashMatchAlgo::GenerateUnfittedFlashHypothesis(constant_qcluster);

	      // At this point I have a TypeError - I will need 'Flash_t' data, but all I have are the 'opflash' data for each track.  There is a function for going from Flash->OpFlash, but not vice versa.
	      
		
	       
	      
	    
	    // First, convert the 3D vector into wire and tick values, which can then be converted into pixel values.  I will still look through the coordinates because looping through the pixels would take longer (1 cm separation between grid points but only 0.3 cm separation between wires on each of the planes).
	    
	    

	    
	    
	    
    

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

      bool contained = false;
      if ( bbox[0][0]>=-10 && bbox[0][0]<=275 && bbox[0][1]>=-10 && bbox[0][1]<=275 )
	contained = true;

      if ( inbbox ) {
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

    }
    std::cout << "Number of Vertex Enclosing BBox: " << num_true_bbox << std::endl;

    tree->Fill();

    //if ( ientry>=3)
    //  break;
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
