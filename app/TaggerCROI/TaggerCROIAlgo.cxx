#include "TaggerCROIAlgo.h"

#include <sstream>
#include <ctime>
#include <stdexcept>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "DataFormat/Pixel2D.h"
#include "DataFormat/ROI.h"

// larlitecv
#include "TaggerTypes/dwall.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "UnipolarHack/UnipolarHackAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgoConfig.h"
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/EndPointFilter.h"
#include "ThruMu/RadialEndpointFilter.h"
#include "ThruMu/PushBoundarySpacePoint.h"
#include "ThruMu/ThruMuTracker.h"
#include "StopMu/StopMuFilterSpacePoints.h"
#include "StopMu/StopMuCluster.h"
#include "StopMu/StopMuFoxTrot.h"
#include "UntaggedClustering/ClusterGroupAlgo.h"
#include "UntaggedClustering/ClusterGroupMatchingAlgo.h"
#include "ContainedROI/TaggerFlashMatchAlgo.h"
#include "ContainedROI/TaggerFlashMatchAlgoV2.h"
#include "ContainedROI/MatchTaggerData2Flash.h"
#include "T3DMerge/T3DCluster.h"
#include "T3DMerge/Track3DRecluster.h"
#include "T3DMerge/T3DPCMerge.h"
#include "T3DMerge/T3D2LarliteTrack.h"
#include "TaggerContourTools/CACAEndPtFilter.h"
#include "MCTruthTools/MCTrack2ImagePath.h"

namespace larlitecv {

  TaggerCROIAlgo::TaggerCROIAlgo( const TaggerCROIAlgoConfig& config )
   : m_config(config) {
    m_time_tracker.resize( kNumStages, 0.0 );
  }

  InputPayload TaggerCROIAlgo::loadInput( DataCoordinator& dataco_input  ) {

    InputPayload input;
    input.clear();

    // Get run, subrun, event, entry
    dataco_input.get_id(input.run,input.subrun,input.event);
    input.entry = dataco_input.get_larcv_io().current_entry();

    // get images (from larcv)
    larcv::EventImage2D* event_imgs = NULL;
    try {
      event_imgs    = (larcv::EventImage2D*)dataco_input.get_larcv_data( larcv::kProductImage2D, m_config.larcv_image_producer );
      if ( event_imgs->Image2DArray().size()==0) {
	if ( !m_config.skip_empty_events )
	  throw std::runtime_error("Number of images=0. LArbys.");
	else {
	  std::cout << "Skipping Empty Events." << std::endl;
	  throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
	}
      }
      else if ( event_imgs->Image2DArray().size()!=3 ) {
	throw std::runtime_error("Number of Images!=3. Weird.");
      }
    }
    catch (const std::exception &exc) {
      std::cerr << "Error retrieving TPC Images: " << exc.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }

    // ------------------------------------------------------------------------------------------//
    // MODIFY IMAGES
    input.img_v.clear();
    event_imgs->Move( input.img_v );
    try {

      if ( m_config.DeJebWires ) {
	for ( auto &img : input.img_v ) {
	  const larcv::ImageMeta& meta = img.meta();
	  for ( int col=0; col<(int)meta.cols(); col++) {
	    for (int row=0; row<(int)meta.rows(); row++) {
	      float val = img.pixel( row, col );
	      if (meta.plane()==0 ) {
		if ( (col>=2016 && col<=2111) || (col>=2176 && 2212) ) {
		  val *= m_config.jebwiresfactor;
		}
	      }
	      img.set_pixel(row,col,val);
	    }
	  }
	}
      }
    }
    catch (const std::exception& e ) {
      std::cerr << "Error Dejebbing wires: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }


    larlitecv::UnipolarHackAlgo unihackalgo;
    try {
      if ( m_config.apply_unipolar_hack ) {
	std::vector<int> applyhack(3,0);
	applyhack[1] = 1;
	std::vector<float> hackthresh(3,-10.0);
	input.img_v = unihackalgo.hackUnipolarArtifact( input.img_v, applyhack, hackthresh );
      }
      if ( input.img_v.size()!=3 )
	throw std::runtime_error("Number of Images incorrect.");
    }
    catch ( const std::exception& e ) {
      std::cerr << "Error making unipolar artifact: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }

    // ------------------------------------------------------------------------------------------//
    // LABEL EMPTY CHANNELS
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::vector< larcv::Image2D > emptyimgs;
    try {
      for ( auto const &img : event_imgs->Image2DArray() ) {
	int p = img.meta().plane();
	larcv::Image2D emptyimg = emptyalgo.labelEmptyChannels( m_config.emptych_thresh.at(p), img );
	emptyimgs.emplace_back( emptyimg );
      }
    }
    catch ( const std::exception& e ) {
      std::cerr << "Error making empty channels: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");      
    }

    // ------------------------------------------------------------------------------------------//
    // LABEL BAD CHANNELS

    input.badch_v.clear();
    try {
      if ( m_config.chstatus_datatype=="LARLITE" ) {
	larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco_input.get_larlite_data( larlite::data::kChStatus, "chstatus" );
	input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
	ev_status->clear(); // clear, we copied the info
      }
      else if ( m_config.chstatus_datatype=="LARCV" ) {
	larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco_input.get_larcv_data( larcv::kProductChStatus, m_config.larcv_chstatus_producer );
	input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
	ev_status->clear(); // clear, we copied the info
      }
      else if ( m_config.chstatus_datatype=="NONE" ) {
	for ( auto const& img : input.img_v ) {
	  larcv::Image2D badch( img.meta() );
	  badch.paint(0.0);
	  input.badch_v.emplace_back( std::move(badch) );
	}
      }
      else {
	throw std::runtime_error("ERROR: ChStatusDataType must be either LARCV or LARLITE");
      }
      
      if ( input.badch_v.size()!=3 ) {
	throw std::runtime_error("Number of Bad Channels not correct.");
      }
    }
    catch ( const std::exception& e ) {
      std::cerr << "Error making badch image: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }

    // ------------------------------------------------------------------------------------------//
    // MAKE GAP CHANNEL IMAGE

    int maxgap = 200; // this seems bad to have this parameter here and not in a config file
    std::vector< larcv::Image2D> gapchimgs_v;
    try {
      gapchimgs_v = emptyalgo.findMissingBadChs( input.img_v, input.badch_v, 5, maxgap );
      // combine with badchs
      for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
        larcv::Image2D& gapchimg = gapchimgs_v[p];
        gapchimg += input.badch_v[p];
      }

      // Set BadCh in input data
      for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
	input.gapch_v.emplace_back( std::move(gapchimgs_v.at(p)) );
      }
      if ( input.gapch_v.size()!=3 )
    	throw std::runtime_error("Number of Gap Channels not correct.");
    }
    catch ( const std::exception& e ) {
      std::cerr << "Error making gapch images: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");      
    }

    // -------------------------------------------------------------------------------------------//
    // COLLECT FLASH INFORMATION
    
    try {
      for ( auto &flashproducer : m_config.opflash_producers ) {
        larlite::event_opflash* opdata = (larlite::event_opflash*)dataco_input.get_larlite_data(larlite::data::kOpFlash, flashproducer );
        std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
        input.opflashes_v.push_back( opdata );
      }
    }
    catch ( const std::exception& e ) {
      std::cerr << "Error retrieving flash information: " << e.what() << std::endl;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }

    // -------------------------------------------------------------------------------------------//
    // LOAD MC TRACK INFORMATION: USED FOR PERFORMANCE METRICS
    if ( m_config.load_mctrack ) {
      try {
	input.p_ev_mctrack = (larlite::event_mctrack*)dataco_input.get_larlite_data(larlite::data::kMCTrack, m_config.mctrack_producer);
      }
      catch (const std::exception& e ) {
	std::cerr << "Error retrieving MC track information upon request: " << e.what() << std::endl;
	throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }
    }

    // -------------------------------------------------------------------------------------------//
    // LOAD EVENT TRIGGER
    try {
      input.p_ev_trigger = (larlite::trigger*)dataco_input.get_larlite_data(larlite::data::kTrigger, m_config.trigger_producer);
    }
    catch (const std::exception& e ) {
      std::cerr << "Error retrieving MC track information upon request: " << e.what() << std::endl;
      input.p_ev_trigger = NULL;
      throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
    }    

    return input;
  }
  
  ThruMuPayload TaggerCROIAlgo::runThruMu( const InputPayload& input ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run ThruMu ==================================" << std::endl;

    ThruMuPayload output;


    if ( !m_config.use_truth_endpoints )
      runBoundaryTagger( input, output );
    else
      runTruthBoundaryTagger( input, output );  // end points using MC truth tracks      

    if ( m_config.use_truth_muontracks ) {
      runTruthThruMu( input, output );
    }
    else {
      runThruMuTracker( input, output);
    }

    if ( m_config.verbosity>0 )
      std::cout << "== End of ThruMu ===============================" << std::endl;

    // return stage output
    return output;
  }

  void TaggerCROIAlgo::runBoundaryTagger( const InputPayload& input, ThruMuPayload& output ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run Boundary Tagger ===============================" << std::endl;

    // configure different stages of the Thrumu Tagger
    std::clock_t timer;

    // (0) make contours
    timer = std::clock();
    m_bmtcv_algo.analyzeImages( input.img_v, input.badch_v, 10.0, 2 );
    m_time_tracker[kThruMuContour] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;
    
    // (1) side tagger
    timer = std::clock();
    larlitecv::BoundaryMuonTaggerAlgo sidetagger;
    sidetagger.configure( m_config.sidetagger_cfg );
    //sidetagger.printConfiguration();

    // (2) flash tagger
    larlitecv::FlashMuonTaggerAlgo anode_flash_tagger(   larlitecv::FlashMuonTaggerAlgo::kAnode );
    larlitecv::FlashMuonTaggerAlgo cathode_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kCathode );
    larlitecv::FlashMuonTaggerAlgo imgends_flash_tagger( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );

    anode_flash_tagger.configure(   m_config.flashtagger_cfg );
    cathode_flash_tagger.configure( m_config.flashtagger_cfg );
    imgends_flash_tagger.configure( m_config.flashtagger_cfg );

    // loading time
    m_time_tracker[kThruMuConfig] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;

    // RUN THE THRUMU ALGOS

    // run side tagger
    timer = std::clock();
    sidetagger.searchforboundarypixels3D( input.img_v, input.badch_v, output.side_spacepoint_v, output.boundarypixel_image_v, output.realspacehit_image_v );
    int nsides[4] = {0};
    for ( auto const& sp : output.side_spacepoint_v ) {
      nsides[ sp.at(0).type ]++;
    }
    if ( m_config.verbosity>=0 ) {
      std::cout << " Side Tagger End Points: " << output.side_spacepoint_v.size() << std::endl;
      std::cout << "   Top: "        << nsides[0] << std::endl;
      std::cout << "   Bottom: "     << nsides[1] << std::endl;
      std::cout << "   Upstream: "   << nsides[2] << std::endl;
      std::cout << "   Downstream: " << nsides[3] << std::endl;
    }
    m_time_tracker[kThruMuBMT] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // run flash tagger
    timer = std::clock();
    anode_flash_tagger.flashMatchTrackEnds(   input.opflashes_v, input.img_v, input.badch_v, output.anode_spacepoint_v );
    cathode_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.cathode_spacepoint_v  );
    imgends_flash_tagger.findImageTrackEnds( input.img_v, input.badch_v, output.imgends_spacepoint_v  );

    int totalflashes = (int)output.anode_spacepoint_v.size() + (int)output.cathode_spacepoint_v.size() + (int)output.imgends_spacepoint_v.size();
    if ( m_config.verbosity>=0 ) {
      std::cout << " Flash Tagger End Points: " << totalflashes << std::endl;
      std::cout << "  Anode: "      << output.anode_spacepoint_v.size() << std::endl;
      std::cout << "  Cathode: "    << output.cathode_spacepoint_v.size() << std::endl;
      std::cout << "  Image Ends: " << output.imgends_spacepoint_v.size() << std::endl;
    }
    std::cout << "Anode spacepoint flash indices: " << std::endl;
    for (int i=0; i<(int)output.anode_spacepoint_v.size(); i++) {
      std::cout << "    [" << i << "] flashidx="
		<< "(" << output.anode_spacepoint_v.at(i).getFlashIndex().ivec << ","
		<< output.anode_spacepoint_v.at(i).getFlashIndex().idx << ")"	
		<< std::endl;
    }
    std::cout << "Cathode spacepoint flash indices: " << std::endl;
    for (int i=0; i<(int)output.cathode_spacepoint_v.size(); i++) {
      std::cout << "    [" << i << "] flashidx="
		<< "(" << output.cathode_spacepoint_v.at(i).getFlashIndex().ivec << ","
		<< output.cathode_spacepoint_v.at(i).getFlashIndex().idx << ")"	
		<< std::endl;
    }
    
    m_time_tracker[kThruMuFlash] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // run end point filters

    //larlitecv::RadialEndpointFilter radialfilter;  // remove end points that cannot form a 3d segment nearby [deprecated]
    //larlitecv::PushBoundarySpacePoint endptpusher; // remove endpt [deprecated]
    //larlitecv::EndPointFilter endptfilter; // removes duplicates [deprecated]
    timer = std::clock();
    larlitecv::CACAEndPtFilter cacaalgo;
    cacaalgo.setVerbosity(0);

    // we collect pointers to all the end points (make a copy for now)
    std::vector< larlitecv::BoundarySpacePoint > all_endpoints;

    // gather endpoints from space points
    for (int isp=0; isp<(int)output.side_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.side_spacepoint_v.at( isp ));
      all_endpoints.push_back( *pts );
    }
    for (int isp=0; isp<(int)output.anode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.anode_spacepoint_v.at(isp));
      all_endpoints.push_back( *pts );
    }
    for (int isp=0; isp<(int)output.cathode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.cathode_spacepoint_v.at(isp));
      all_endpoints.push_back( *pts );
    }
    for (int isp=0; isp<(int)output.imgends_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.imgends_spacepoint_v.at(isp));
      all_endpoints.push_back( *pts );
    }
    if ( m_config.verbosity>0 )
      std::cout << "number of endpoints pre-filters: " << all_endpoints.size() << std::endl;

    std::vector< const std::vector<larlitecv::BoundarySpacePoint>* > sp_v;
    sp_v.push_back( &all_endpoints );
    std::vector< std::vector<int> > caca_results;    
    cacaalgo.evaluateEndPoints( sp_v, input.opflashes_v, input.img_v, input.badch_v, m_bmtcv_algo.m_plane_atomicmeta_v, 150.0, caca_results );

    // prepare the boundary points that pass
    std::vector<larlitecv::BoundarySpacePoint> cacapassing_moved_v = cacaalgo.regenerateFitleredBoundaryPoints( input.img_v );
    
    // clean up
    all_endpoints.clear();
    sp_v.clear();
    m_time_tracker[kThruMuFilter] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;
    

    // collect output
    // remove the filtered end points
    for ( size_t idx=0; idx<cacapassing_moved_v.size(); idx++ ) {
      larlitecv::BoundarySpacePoint& sp = cacapassing_moved_v[idx];

      if (sp.type()<=larlitecv::kDownstream ) {
	output.side_filtered_v.emplace_back( std::move(sp) );
      }
      else if (sp.type()==larlitecv::kAnode) {
	output.anode_filtered_v.emplace_back( std::move(sp) );
      }
      else if (sp.type()==larlitecv::kCathode) {
	output.cathode_filtered_v.emplace_back( std::move(sp) );
      }
      else if (sp.type()==larlitecv::kImageEnd) {
	output.imgends_filtered_v.emplace_back( std::move(sp) );
      }
      else {
	std::stringstream ss;
	ss << __FILE__ << ":" << __LINE__ << " unrecognized boundary type" << std::endl;
	throw std::runtime_error(ss.str());
      }
    }

    if ( m_config.verbosity>=0 ) {
      std::cout << " Filtered Side Tagger End Points: " << cacapassing_moved_v.size() << std::endl;
      std::cout << "   Side: "        << output.side_filtered_v.size() << std::endl;
      std::cout << "   Anode: "       << output.anode_filtered_v.size() << std::endl;
      std::cout << "   Cathode: "     << output.cathode_filtered_v.size() << std::endl;
      std::cout << "   ImageEnds: "   << output.imgends_filtered_v.size() << std::endl;
    }

    
    if ( m_config.verbosity>0 )
      std::cout << "== End of Boundary Tagger ===============================" << std::endl;

    return;
  }

  void TaggerCROIAlgo::runTruthBoundaryTagger( const InputPayload& input, ThruMuPayload& output ) {
    // This uses MC truth information to make boundary end points
    // We use this for debugging and performance metrics

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    const larcv::ImageMeta& meta = input.img_v.front().meta();

    m_truthxingdata.clear();
    larlitecv::analyzeCrossingMCTracks( m_truthxingdata, meta, input.img_v, input.p_ev_trigger, input.p_ev_mctrack, input.opflashes_v, false );

    // extract the truth points and make end point containers
    for ( auto const& xingpt : m_truthxingdata.truthcrossingptinfo_v ) {
      std::vector<float> tyz(3);
      for (int i=0; i<3; i++)
	tyz[i] = xingpt.crossingpt_detsce_tyz[i];
      tyz[0] = (tyz[0]-3200.0)*cm_per_tick;
      BoundarySpacePoint sp( (larlitecv::BoundaryEnd_t)xingpt.type, tyz, meta );
      if ( sp.type()<=larlitecv::kDownstream  ) {
	output.side_spacepoint_v.push_back( sp );
	output.side_filtered_v.push_back( sp );	
      }
      else if ( sp.type()==larlitecv::kAnode ) {
	output.anode_spacepoint_v.push_back( sp );
	output.anode_filtered_v.push_back( sp );		
      }
      else if ( sp.type()==larlitecv::kCathode ) {
	output.cathode_spacepoint_v.push_back( sp );
	output.cathode_filtered_v.push_back( sp );
      }
      else if ( sp.type()==larlitecv::kImageEnd ) {
	output.imgends_spacepoint_v.push_back( sp );
	output.imgends_filtered_v.push_back( sp );	
      }
    }
    
  }

  void TaggerCROIAlgo::runThruMuTracker( const InputPayload& input, ThruMuPayload& output ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run ThruMu Tracker ===============================" << std::endl;

    // configure different stages of the Thrumu Tagger
    std::clock_t timer = std::clock();

    // thrumu tracker
    larlitecv::ThruMuTracker thrumu_tracker( m_config.thrumu_tracker_cfg );

    m_time_tracker[kThruMuConfig] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;

    // RUN THE THRUMU ALGOS

    // run side tagger
    timer = std::clock();

    // we collect pointers to all the end points
    std::vector< const larlitecv::BoundarySpacePoint* > all_endpoints;

    // Collect the indices for the flash that determines each of the endpoints and the index for the producer that makes each of the endpoints as well.
    // This is a precaution to ensure that we are collecting this information properly.
    std::vector< int > all_endpoints_flash_idx_v;
    all_endpoints_flash_idx_v.clear();

    std::vector< int > all_endpoints_boundary_type_idx_v;
    all_endpoints_boundary_type_idx_v.clear();

    // When filling the 'all_endpoints_flash_idx_v' and 'all_endpoints_boundary_type_idx_v',
    // the correct flash indices have already been matched to the correct endpoint.
    // The important thing is that these two vectors are filled at the same point as their corresponding endpoint in 'all_endpoints'.

    // gather endpoints from space points
    for (int isp=0; isp<(int)output.side_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.side_filtered_v.at( isp ));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.anode_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.anode_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.cathode_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.cathode_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.imgends_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.imgends_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    if ( m_config.verbosity>0 )
      std::cout << "number of endpoints to search for thrumu: " << all_endpoints.size() << std::endl;

    // make track clusters
    std::vector<int> used_endpoints( all_endpoints.size(), 0 );
    if ( m_config.run_thrumu_tracker ) {
      timer = std::clock();
      output.trackcluster3d_v.clear();
      output.tagged_v.clear();
      
      thrumu_tracker.makeTrackClusters3D( m_config.croi_selection_cfg.genflashmatch_cfg, input.img_v, input.gapch_v,
					  all_endpoints, output.trackcluster3d_v, 
      					  output.tagged_v, used_endpoints, input.opflashes_v ); 

      m_time_tracker[kThruMuTracker]  +=  (std::clock()-timer)/(double)CLOCKS_PER_SEC;
      if ( m_config.verbosity>0 )
        std::cout << "thrumu tracker search " << all_endpoints.size() << " end points in " << m_time_tracker[kThruMuTracker] << " sec" << std::endl;
    }
    else {
      if ( m_config.verbosity>0 )
        std::cout << "config tells us to skip thrumu track." << std::endl;
    }
    
    // collect unused endpoints
    output.used_spacepoint_v.clear();
    output.unused_spacepoint_v.clear();
    for ( size_t isp=0; isp<all_endpoints.size(); isp++ ) {
      if ( used_endpoints.at(isp)==1 )
        output.used_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
      else {
        const BoundarySpacePoint& sp = *(all_endpoints.at(isp));
        if ( m_config.verbosity>1 )
          std::cout << "unused spacepoint for StopMu: (" << sp.pos()[0] << "," << sp.pos()[1] << "," << sp.pos()[2] << ")" << std::endl;
        output.unused_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
      }
    }

    // copy track and pixels into separate containers.
    output.track_v.clear();
    output.pixelcluster_v.clear();
    for ( auto const& bmtrack : output.trackcluster3d_v ) {
      output.track_v.push_back( bmtrack.makeTrack() );
      std::vector< larcv::Pixel2DCluster > cluster_v;
      for ( auto const& track2d : bmtrack.plane_pixels )
        cluster_v.push_back( track2d );
      output.pixelcluster_v.emplace_back( std::move(cluster_v) );
    }

    if ( m_config.verbosity>0 )
      std::cout << "== End of ThruMu Tracker ===============================" << std::endl;

  }

  void TaggerCROIAlgo::runTruthThruMu( const InputPayload& input, ThruMuPayload& output ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run MC Truth ThruMu ===============================" << std::endl;

    // We fill the ThruMu payload using MC Track Truth

    std::clock_t timer = std::clock();

    // we collect pointers to all the end points
    std::vector< const larlitecv::BoundarySpacePoint* > all_endpoints;

    // Collect the indices for the flash that determines each of the endpoints and the index for the producer that makes each of the endpoints as well.
    // This is a precaution to ensure that we are collecting this information properly.
    std::vector< int > all_endpoints_flash_idx_v;
    all_endpoints_flash_idx_v.clear();

    std::vector< int > all_endpoints_boundary_type_idx_v;
    all_endpoints_boundary_type_idx_v.clear();

    // gather endpoints from space points
    for (int isp=0; isp<(int)output.side_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.side_filtered_v.at( isp ));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.anode_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.anode_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.cathode_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.cathode_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)output.imgends_filtered_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(output.imgends_filtered_v.at(isp));
      all_endpoints.push_back( pts );
    }
    if ( m_config.verbosity>0 )
      std::cout << "number of endpoints to search for thrumu: " << all_endpoints.size() << std::endl;

    // make track clusters
    std::vector<int> used_endpoints( all_endpoints.size(), 1 ); // all end points "used"
    output.trackcluster3d_v.clear();
    output.tagged_v.clear();
    
    // we loop through the true crossing points saved in m_truthxingdata
    // we find those tagged as through-going
    for ( auto const& truthxing : m_truthxingdata.truthcrossingptinfo_v ) {
      if ( truthxing.start_or_end==1 )
	continue; // we only want to work with starts

      // get the mc track
      int mctrackidx = truthxing.mctrack_index;
      const larlite::mctrack& truthtrack = input.p_ev_mctrack->at( mctrackidx );

      // get the other end point
      const TruthCrossingPointAna_t* endpt = NULL;
      for ( auto const& truthendpt : m_truthxingdata.truthcrossingptinfo_v ) {
	if ( truthendpt.start_or_end==0 )
	  continue;
	if ( truthendpt.mctrack_index==mctrackidx ) {
	  endpt = &truthendpt;
	  break;
	}
      }
      
      // use the above to build the BMTrackCluster3D object!!
      // make path: we have to go between start and end pt of mctrack object (just use the whole thing for, restricting inside TPC)
      std::vector<std::vector<double>> path_tyz_sce = mctrack2tyz( truthtrack, input.p_ev_trigger->TriggerTime(), true, &m_sce );
      
      BoundarySpacePoint startspt( (BoundaryEnd_t)truthxing.type, truthxing.crossingpt_detsce, input.img_v.front().meta() );
      BoundarySpacePoint endspt;

      if ( endpt!=NULL )
	endspt = BoundarySpacePoint( (BoundaryEnd_t)endpt->type, endpt->crossingpt_detsce, input.img_v.front().meta() );

      BMTrackCluster3D track3d( startspt, endspt, path_tyz_sce );

      /// make pixels and tag
      track3d.markImageWithTrack( input.img_v, input.badch_v,
				  m_config.thrumu_tracker_cfg.pixel_threshold,
				  m_config.thrumu_tracker_cfg.tag_neighborhood,
				  output.tagged_v, 0.3, 10.0 );
      track3d.mctrackid = truthxing.mctrackid;
      std::cout << "Fill MCTrack trackindex=" << truthxing.mctrack_index << " mctrackid=" << track3d.mctrackid << std::endl;
      output.trackcluster3d_v.emplace_back( std::move(track3d) );
    }
    
    // collect unused endpoints
    output.used_spacepoint_v.clear();
    output.unused_spacepoint_v.clear();
    for ( size_t isp=0; isp<all_endpoints.size(); isp++ ) {
      if ( used_endpoints.at(isp)==1 )
        output.used_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
      else {
        const BoundarySpacePoint& sp = *(all_endpoints.at(isp));
        if ( m_config.verbosity>1 )
          std::cout << "unused spacepoint for StopMu: (" << sp.pos()[0] << "," << sp.pos()[1] << "," << sp.pos()[2] << ")" << std::endl;
        output.unused_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
      }
    }

    // copy track and pixels into separate containers.
    output.track_v.clear();
    output.pixelcluster_v.clear();
    for ( auto const& bmtrack : output.trackcluster3d_v ) {
      output.track_v.push_back( bmtrack.makeTrack() );
      std::vector< larcv::Pixel2DCluster > cluster_v;
      for ( auto const& track2d : bmtrack.plane_pixels )
        cluster_v.push_back( track2d );
      output.pixelcluster_v.emplace_back( std::move(cluster_v) );
    }
    
    if ( m_config.verbosity>0 )
      std::cout << "== End of Truth ThruMu ===============================" << std::endl;

    m_time_tracker[kThruMuTracker]  +=  (std::clock()-timer)/(double)CLOCKS_PER_SEC;
  }
  
  StopMuPayload TaggerCROIAlgo::runStopMu( const InputPayload& input, const ThruMuPayload& thrumu ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run StopMu Tracker ===============================" << std::endl;

    StopMuPayload output;

    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;

    // Algos
    larlitecv::StopMuFilterSpacePoints stopmu_filterpts( m_config.stopmu_filterpts_cfg );
    larlitecv::StopMuCluster           stopmu_cluster( m_config.stopmu_cluster_cfg );
    larlitecv::StopMuFoxTrot           stopmu_foxtrot( m_config.stopmu_foxtrot_cfg );

    if ( m_config.stopmu_cluster_cfg.save_pass_images || m_config.stopmu_cluster_cfg.dump_tagged_images  ) {
      std::stringstream ss;
      ss << "smc_r" << input.run << "_s" << input.subrun << "_e" << input.event << "_i" << input.entry;
      stopmu_cluster.setOpenCVImageStemName( ss.str() );
    }

    // We strip the unused pixel locations into vector of pixels.
    std::vector< std::vector< const larcv::Pixel2D* > > unused_spacepoint_v;

    const larcv::ImageMeta& meta = input.img_v.front().meta();

    for ( auto const& pt : thrumu.unused_spacepoint_v ) {

      Double_t xyz[3];
      for (int i=0; i<3; i++)
        xyz[i] = pt.pos()[i];

      float tick = xyz[0]/cm_per_tick + 3200.0;
      if ( tick<=meta.min_y() ) tick = meta.min_y()+1;
      if ( tick>=meta.max_y() ) tick = meta.max_y()-1;
      int row = meta.row( tick );
      std::vector<float> wid( input.img_v.size(),0.0);
      for (size_t p=0; p<input.img_v.size(); p++) {
        float wire = larutil::Geometry::GetME()->WireCoordinate(xyz,p);
        if ( wire<0 ) wire = 0;
        if ( wire>=3456 ) wire = 3455;
        int col = input.img_v.at(p).meta().col( wire );
        larcv::Pixel2D pix( col, row );
        output.stopmu_pixel_endpt_v.Emplace( p, std::move(pix) );
      }

    }

    std::vector< larcv::EventPixel2D* > unused_spacepoints;
    unused_spacepoints.push_back( &(output.stopmu_pixel_endpt_v) );

    output.stopmu_candidate_endpt_v = stopmu_filterpts.filterSpacePoints( unused_spacepoints, thrumu.tagged_v, input.badch_v );
    if ( m_config.verbosity>0 )
      std::cout << "  Number of candidate stop-mu start points: " << output.stopmu_candidate_endpt_v.size() << std::endl;

    std::clock_t timer = std::clock();
    //output.stopmu_trackcluster_v = stopmu_cluster.findStopMuTracks( input.img_v, input.gapch_v, thrumu.tagged_v, output.stopmu_candidate_endpt_v );
    output.stopmu_trackcluster_v = stopmu_foxtrot.findStopMuTracks( input.img_v, input.gapch_v, thrumu.tagged_v, thrumu.unused_spacepoint_v );
    m_time_tracker[kStopMuTracker] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;
    if ( m_config.verbosity>0 )
      std::cout << "  Number of candidate StopMu tracks: " << output.stopmu_trackcluster_v.size() << std::endl;
    
    // make stopmu-tagged image
    for (size_t p=0; p<input.img_v.size(); p++) {
      larcv::Image2D stopmu_img( input.img_v.at(p).meta() );
      stopmu_img.paint(0);
      output.stopmu_v.emplace_back( std::move(stopmu_img) );
    }
    for ( size_t itrack=0; itrack<output.stopmu_trackcluster_v.size(); itrack++ ) {
      larlitecv::BMTrackCluster3D& track3d = output.stopmu_trackcluster_v.at(itrack);
      std::vector<larcv::Pixel2DCluster> trackpix_v = track3d.getTrackPixelsFromImages( input.img_v, input.gapch_v,
											m_config.stopmu_foxtrot_cfg.foxtrotalgo_cfg.pixel_thresholds,
											m_config.thrumu_tracker_cfg.tag_neighborhood, 0.3 );
      for (size_t p=0; p<trackpix_v.size(); p++) {
        const larcv::Pixel2DCluster& trackpixs = trackpix_v[p];
        for ( auto const& pix : trackpixs ) {
          output.stopmu_v.at(p).set_pixel( pix.Y(), pix.X(), 255 );
        }
      }
    }
    
    // copy track and pixels into separate containers.
    for ( auto const& bmtrack : output.stopmu_trackcluster_v ) {
      output.track_v.push_back( bmtrack.makeTrack() );
      std::vector< larcv::Pixel2DCluster > cluster_v;
      for ( auto const& track2d : bmtrack.plane_pixels )
        cluster_v.push_back( track2d );
      output.pixelcluster_v.emplace_back( std::move(cluster_v) );
    }
    
    if ( m_config.verbosity>0 )
      std::cout << "== End of StopMu Tracker ===============================" << std::endl;
    
    return output;
  }

  CROIPayload TaggerCROIAlgo::runCROISelection( const InputPayload& input, const ThruMuPayload& thrumu, const StopMuPayload& stopmu ) {

    if ( m_config.verbosity>0 )
      std::cout << "== Run Untagged Cluster and CROI Selection ===============================" << std::endl;

    CROIPayload output;

    ClusterGroupAlgo         clusteralgo(   m_config.untagged_cluster_cfg );
    ClusterGroupMatchingAlgo matchingalgo;
    TaggerFlashMatchAlgo     selectionalgo(   m_config.croi_selection_cfg );
    TaggerFlashMatchAlgoV2   selectionv2algo( m_config.croi_selection_cfg );
    Track3DRecluster         reclusteralgo;
    T3DPCMerge               pcamergealgo;
    //pcamergealgo.setVerbosity(0);

    std::clock_t timer;

    if ( m_config.recluster_stop_and_thru ) {

      timer = std::clock();

      // we recluster tracks
      for ( auto const& llthrumu : thrumu.track_v ) {
        int npts = llthrumu.NumberTrajectoryPoints();
        std::vector< std::vector<float> > path;
        for (int n=0; n<npts; n++) {
          std::vector<float> pos(3);
          for (int v=0; v<3; v++)
            pos[v] = llthrumu.LocationAtPoint(n)[v];
          path.emplace_back( std::move(pos) );
        }
        reclusteralgo.addPath( path );
      }
      for ( auto const& llstopmu : stopmu.track_v ) {
        int npts = llstopmu.NumberTrajectoryPoints();
        std::vector< std::vector<float> > path;
        for (int n=0; n<npts; n++) {
          std::vector<float> pos(3);
          for (int v=0; v<3; v++)
          pos[v] = llstopmu.LocationAtPoint(n)[v];
          path.emplace_back( std::move(pos) );
        }
        reclusteralgo.addPath( path );
      }

      if ( m_config.verbosity>0 )
        std::cout << "Run Recluster Algo." << std::endl;
      output.stopthru_reclustered_v = reclusteralgo.recluster();

      // Use reclustering to tag image
      for ( size_t p=0; p<input.img_v.size(); p++) {
        larcv::Image2D tagged( input.img_v.at(p).meta() );
        tagged.paint(0.0);
        larcv::Image2D sub( input.img_v.at(p) );
        output.tagged_v.emplace_back( std::move(tagged) );
        output.subimg_v.emplace_back( std::move(sub) );
      }

      for ( auto& t3dtrack : output.stopthru_reclustered_v ) {

        std::vector< larcv::Pixel2DCluster > pixel_v = t3dtrack.getPixelsFromImages( input.img_v, input.gapch_v,
        									     m_config.thrumu_tracker_cfg.pixel_threshold,
        									     m_config.thrumu_tracker_cfg.tag_neighborhood, 0.3 );
        for (size_t p=0; p<pixel_v.size(); p++) {
          larcv::Image2D& tagged = output.tagged_v[p];
          larcv::Image2D& sub    = output.subimg_v[p];
          for ( auto const& pix : pixel_v[p] ) {
            tagged.set_pixel(pix.Y(),pix.X(),255);
            // subtraction image: below threshold and tagged pixels get zeroed (for clustering)
            if ( sub.pixel(pix.Y(),pix.X())<m_config.thrumu_tracker_cfg.pixel_threshold[p] )
              sub.set_pixel(pix.Y(),pix.X(),0.0);
          }
        }
        output.stopthru_reclustered_pixels_v.emplace_back( std::move(pixel_v) );
      }// loop over treclustered tracks

      m_time_tracker[kRecluster] = (std::clock()-timer)/(double)CLOCKS_PER_SEC;
    }
    else {

      m_time_tracker[kRecluster] = 0.0; // not run

      // Make Thrumu/StopMu tagged and subtracted images
      //std::vector< larcv::Image2D > tagged_v;
      //std::vector< larcv::Image2D > subimg_v;
      for ( size_t p=0; p<input.img_v.size(); p++) {
        larcv::Image2D tagged( input.img_v.at(p).meta() );
        tagged.paint(0.0);
        larcv::Image2D sub( input.img_v.at(p) );

        for ( size_t r=0; r<tagged.meta().rows(); r++ ) {
          for ( size_t c=0; c<tagged.meta().cols(); c++ ) {
            // tagged image
            if ( thrumu.tagged_v.at(p).pixel(r,c)>0 || stopmu.stopmu_v.at(p).pixel(r,c)>0 )
              tagged.set_pixel(r,c,255);
            // subtraction image: below threshold and tagged pixels get zeroed (for clustering)
            if ( sub.pixel(r,c)<10.0 || thrumu.tagged_v.at(p).pixel(r,c)>0 || stopmu.stopmu_v.at(p).pixel(r,c)>0 )
              sub.set_pixel(r,c,0.0);
          }
        }
        output.tagged_v.emplace_back( std::move(tagged) );
        output.subimg_v.emplace_back( std::move(sub) );
      }
    }//end of prep when NOT reclustering stop/mu

    // ----------------------
    //  RUN ALGOS

    timer = std::clock();

    output.plane_groups_v = clusteralgo.MakeClusterGroups( input.img_v, input.gapch_v, output.tagged_v );
    //output.plane_groups_v = clusteralgo.MakeClusterGroups( output.subimg_v, input.gapch_v, output.tagged_v );

    output.vols_v = matchingalgo.MatchClusterGroups( output.subimg_v, output.plane_groups_v );

    m_time_tracker[kUntagged] = (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // ------------------------------------------------------------------
    // WE COLLECT OUR CLUSTER DATA, forming TaggerFlashMatchData objects
    // Collect Pixel2D clusters for each plane. Define LArLite Track
    //std::vector< larlitecv::TaggerFlashMatchData > flashdata_v;

    if ( m_config.recluster_stop_and_thru ) {
      // USE RECLUSTERED STOP/THRU TRACKS
      /*
      for (int itrack=0; itrack<(int)output.stopthru_reclustered_v.size(); itrack++) {
        larlite::track lltrack = larlitecv::T3D2LarliteTrack( output.stopthru_reclustered_v[itrack] );
        larlitecv::TaggerFlashMatchData reclustered_track( larlitecv::TaggerFlashMatchData::kThruMu, output.stopthru_reclustered_pixels_v[itrack], lltrack );
        output.flashdata_v.emplace_back( std::move(reclustered_track) );
      }
      // ASSOCIATE FLASHES TO THE TAGGER DATA
      larlitecv::MatchTaggerData2Flash( output.flashdata_v, input.opflashes_v, thrumu.anode_spacepoint_v, thrumu.cathode_spacepoint_v, 10.0 );
      */
    }
    else {
      // DONT USE RECLUSTERED STOP/THRU MU TRACKS
      // ThruMu
      for ( int itrack=0; itrack<(int)thrumu.trackcluster3d_v.size(); itrack++ ) {
	
        larlitecv::TaggerFlashMatchData thrumu_track( larlitecv::TaggerFlashMatchData::kThruMu, thrumu.pixelcluster_v.at(itrack), thrumu.track_v.at(itrack) );
	thrumu_track.mctrackid = thrumu.trackcluster3d_v[itrack].mctrackid;
	std::cout << "transfer mctrackid=" << thrumu_track.mctrackid << " " << thrumu.trackcluster3d_v[itrack].mctrackid << std::endl;
        output.flashdata_v.emplace_back( std::move(thrumu_track) );
      }

      // StopMu
      for ( int itrack=0; itrack<(int)stopmu.stopmu_trackcluster_v.size(); itrack++ ) {
        larlitecv::TaggerFlashMatchData stopmu_track( larlitecv::TaggerFlashMatchData::kStopMu, stopmu.pixelcluster_v.at(itrack), stopmu.track_v.at(itrack) );
	stopmu_track.mctrackid = stopmu.stopmu_trackcluster_v[itrack].mctrackid;	
        output.flashdata_v.emplace_back( std::move(stopmu_track) );
      }
    }

    // ------------------------------------------------------------------
    // FIND CONTAINED CLUSTERS
    timer = std::clock();

    // Find Contained Clusters
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    std::vector<larlitecv::TaggerFlashMatchData> contained_tracks_v;
    for ( auto const& vol : output.vols_v ) {

      // there should be a selection here...
      // select by (1) good fraction (2) charge even-ness
      if ( vol.frac_good_slices<0.8 )
        continue;

      std::cout << "VOL: clgroup[" << vol._clustergroup_indices[0] << "," << vol._clustergroup_indices[1] << "," << vol._clustergroup_indices[2] << "] "
                << " numslices=" << vol.num_slices
                << " goodslices=" << vol.num_good_slices
                << " fracgood=" << vol.frac_good_slices
                << " planecharge=[" << vol.plane_charge[0] << "," << vol.plane_charge[1] << "," << vol.plane_charge[2] << "]"
                << std::endl;

      // we need to make a larlite::track object for this. we use the centroid of the slices
      std::vector< std::vector<float> > xyz_v;
      std::vector< float > slice_charge;
      for ( size_t islice=0; islice<vol.slices.size(); islice++ ) {

        const larlitecv::Slice_t& slice = vol.slices.at(islice);

        if ( slice.inside_tpc_boundary.size()==0 )
          continue;

        // for each slice volume, we are going to do the easy thing first and represent charge at centroid
        // of volume. not great, I know.
        float centroid[2] = {0.0};
        for ( auto const pt : slice.inside_tpc_boundary ) {
          for (int i=0; i<2; i++)
            centroid[i] += pt[i];
        }
        for (int i=0; i<2; i++) {
          centroid[i] /= float(slice.inside_tpc_boundary.size());
        }

        float tick = input.img_v.front().meta().pos_y( 0.5*(slice.row_interval[0]+slice.row_interval[1]) );
        float x = (tick-3200.0)*cm_per_tick;
        std::vector< float > xyz(3,0);
        xyz[0] = x;
        xyz[1] = centroid[0];
        xyz[2] = centroid[1];

        xyz_v.emplace_back( std::move(xyz) );

        // slice charge
        float ave_charge = 0.;
        for ( int i=0; i<3; i++ ) {
          ave_charge += slice.plane_charge[i];
        }
        ave_charge /= 3.0;
        slice_charge.push_back( ave_charge );
      }

      larlite::track contained_track;
      for ( int ipt=0; ipt<(int)xyz_v.size()-1; ipt++ ) {
        const std::vector<float>& xyz = xyz_v.at(ipt);
        TVector3 pos( xyz[0], xyz[1], xyz[2] );
        const std::vector<float>& xyz_next = xyz_v.at(ipt+1);
        float dir[3] = {0};
        float norm = 0.;
        for (int i=0; i<3; i++) {
          dir[i] = xyz_next[i] - xyz[i];
          norm += dir[i]*dir[i];
        }
        norm = sqrt(norm);
        for (int i=0; i<3; i++)
          dir[i] /= norm;
        TVector3 dirv( dir[0], dir[1], dir[2] );
        contained_track.add_vertex( pos );
        contained_track.add_direction( dirv );
        contained_track.add_momentum( slice_charge.at(ipt) );

        if ( ipt==(int)xyz_v.size()-2) {
          TVector3 nextpos( xyz_next[0], xyz_next[1], xyz_next[2] );
          contained_track.add_vertex( nextpos );
          contained_track.add_direction( dirv );
          contained_track.add_momentum( slice_charge.at(ipt+1) );
        }
      }
      if ( contained_track.NumberTrajectoryPoints()==0 )
        continue;
      larlitecv::TaggerFlashMatchData contained_cluster( larlitecv::TaggerFlashMatchData::kUntagged, vol.m_plane_pixels, contained_track );
      //output.flashdata_v.emplace_back( std::move(contained_cluster) );
      contained_tracks_v.emplace_back( std::move(contained_cluster) );
    }// end of vol loop


    m_time_tracker[kUntagged] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // --------------------------------------------------------------------//
    // RECLUSTER CONTAINED TRACKS

    timer = std::clock();

    reclusteralgo.clear();
    for ( auto const& contained : contained_tracks_v ) {
      // convert larlite into T3DCluster (avoiding this copy can come later)
      std::vector< std::vector<double> > path;
      for (int ipt=0; ipt<(int)contained.m_track3d.NumberTrajectoryPoints(); ipt++) {
        std::vector<double> pos(3);
        for (int i=0; i<3; i++)
          pos[i] = contained.m_track3d.LocationAtPoint(ipt)[i];
        path.push_back( pos );
      }
      reclusteralgo.addPath( path );
    }
    std::vector< T3DCluster > reclustered_contained = reclusteralgo.recluster();
    std::cout << "Reclustered Contained Tracks from " << contained_tracks_v.size() << " to " << reclustered_contained.size();

    // Add the contained tracks to the output.
    // Append the contained tracks onto the end of the output.
    output.contained_tracks_v.clear();

    // Place these into the output for the CROI stage.
    for ( auto& contained_track_loop : contained_tracks_v ) {

      output.contained_tracks_v.emplace_back( std::move(contained_track_loop) );

    }

    m_time_tracker[kRecluster] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // --------------------------------------------------------------------//
    // RECLUSTER CONTAINED + STOP/THRUMU
    reclusteralgo.clear();

    for ( auto& t3d : output.stopthru_reclustered_v ) {
      reclusteralgo.addTrack( t3d );
    }
    for ( auto& t3d : reclustered_contained) {
      reclusteralgo.addTrack( t3d );
    }
    std::vector< T3DCluster > reclustered_alltracks_v = reclusteralgo.recluster();

    // Add the reclustered contained tracks to the output.
    // Add the contained reclustered tracks to the output.
    output.reclustered_contained_v.clear();

    // Place these into the output.
    for ( auto& reclustered_contained_track : reclustered_contained ) {

      output.reclustered_contained_v.emplace_back( std::move(reclustered_contained_track) );

    }

    // --------------------------------------------------------------------//
    // PCAMERGE

    timer = std::clock();

    std::vector< T3DCluster > pcmerged_v = pcamergealgo.merge( reclustered_alltracks_v );

    m_time_tracker[kPCAmerge] = (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // --------------------------------------------------------------------//
    // CONVERT TRACKS TO TAGGERFLASHMATCH DATA AND REASSOCIATE THEM TO FLASHES

    for (int itrack=0; itrack<(int)pcmerged_v.size(); itrack++) {
      T3DCluster& t3d = pcmerged_v[itrack];
      if ( t3d.getPath().size()==0 )
        continue;
      larlite::track lltrack = larlitecv::T3D2LarliteTrack( t3d );
      std::vector< larcv::Pixel2DCluster > pixel_v = t3d.getPixelsFromImages( input.img_v, input.gapch_v,
									      m_config.thrumu_tracker_cfg.pixel_threshold,
									      m_config.thrumu_tracker_cfg.tag_neighborhood, 0.3 );
      larlitecv::TaggerFlashMatchData pcmerged_track( larlitecv::TaggerFlashMatchData::kThruMu, pixel_v, lltrack );

      // we tag the type based on dwall
      int start_boundary = -1;
      int end_boundary = -1;
      double start_dwall = larlitecv::dwall( t3d.getPath().front(), start_boundary );
      double end_dwall   = larlitecv::dwall( t3d.getPath().back(), end_boundary );
      if ( start_dwall<10.0 && end_dwall<10.0 ) {
        pcmerged_track.m_type = larlitecv::TaggerFlashMatchData::kThruMu;
      }
      else if ( start_dwall<10.0 || end_dwall<10.0 ) {
        pcmerged_track.m_type = larlitecv::TaggerFlashMatchData::kStopMu;
      }
      else {
        pcmerged_track.m_type = larlitecv::TaggerFlashMatchData::kUntagged;
      }

      output.flashdata_v.emplace_back( std::move(pcmerged_track) );
    }

    // reassociate flashes
    larlitecv::MatchTaggerData2Flash( output.flashdata_v, input.opflashes_v, thrumu.anode_spacepoint_v, thrumu.cathode_spacepoint_v, 10.0 );

    // --------------------------------------------------------------------//
    // SELECT ROIs

    timer = std::clock();

    output.flashdata_selected_v.resize( output.flashdata_v.size(), 0 );
    std::vector<larcv::ROI> selected_rois;
    if ( m_config.croi_selection_cfg.use_version==1 )
      selected_rois = selectionalgo.FindFlashMatchedContainedROIs( output.flashdata_v, input.opflashes_v, output.flashdata_selected_v );
    else if ( m_config.croi_selection_cfg.use_version==2 ) {
      if ( m_config.use_truth_endpoints )
	selectionv2algo.provideTruthCrossingData( m_truthxingdata );
      selected_rois = selectionv2algo.FindFlashMatchedContainedROIs( output.flashdata_v, input.opflashes_v, input.img_v, output.flashdata_selected_v );
    }
    else
      throw std::runtime_error("Unrecognized TaggerFlashMatchAlgo version");
    
    output.croi_v.clear();
    // Moved to croi selection algo
    // for ( size_t itrack=0; itrack<output.flashdata_v.size(); itrack++ ){
    //   const larlite::track& track3d = output.flashdata_v.at(itrack).m_track3d;
    //   if ( output.flashdata_selected_v.at(itrack)==0 || track3d.NumberTrajectoryPoints()==0)
    //     continue;

    //   larcv::ROI croi = output.flashdata_v.at(itrack).MakeROI( input.img_v, m_config.croi_selection_cfg.bbox_pad , true );

    //   std::cout << "[Selected CROI]" << std::endl;
    //   for ( size_t p=0; p<3; p++ ) {
    //     std::cout << "  " << croi.BB(p).dump() << std::endl;
    //   }

    //   output.croi_v.emplace_back( std::move(croi) );
    // }
    for ( auto& roi : selected_rois ) {
      output.croi_v.push_back( roi );
    }

    m_time_tracker[kCROI] = (std::clock()-timer)/(double)CLOCKS_PER_SEC;

    // -----------------------------------------------------------------------
    // Copy cut results
    if ( m_config.croi_write_cfg.get<bool>("WriteCutResults") ) {

      if ( m_config.croi_selection_cfg.use_version==1 ) {
	// version 1
	output.containment_dwall_v = selectionalgo.getContainmentCutValues();
	output.min_chi2_v          = selectionalgo.getInTimeChi2Values();
	output.totpe_peratio_v     = selectionalgo.getTotPEratioValues();
	output.cosmicflash_ratio_dchi_v = selectionalgo.getCosmicDeltaChi2Values();
	
	output.flashdata_passes_containment_v = selectionalgo.getContainmentCutResults();
	output.flashdata_passes_cosmicflash_ratio_v = selectionalgo.getCosmicRatioCutResults();
	output.flashdata_passes_flashmatch_v = selectionalgo.getFlashMatchCutResults();
	output.flashdata_passes_totpe_v      = selectionalgo.getTotalPECutResults();
      }
      else {
	// version 2 (we reuse vectors named for v1 variables)
	output.flashdata_passes_containment_v       = selectionv2algo.getContainmentCutResults();
	output.flashdata_passes_cosmicflash_ratio_v = selectionv2algo.getCosmicFlashMatchResults();
	output.v2_flashdata_passes_flashpos_v       = selectionv2algo.getIntimeFlashPosResults();
	
	output.containment_dwall_v        = selectionv2algo.getContainmentCutValues();
	output.cosmicflash_ratio_dchi_v   = selectionv2algo.getCosmicFlashMatchChi2Values();
	output.cosmicflash_bestindex_v    = selectionv2algo.getCosmicMatchedFlashIndices();
	output.cosmicflash_mctrackid_v    = selectionv2algo.getCosmicMCTrackIDs();
	output.cosmicflash_bestflash_mctrackid_v = selectionv2algo.getCosmicBestFlashMCTrackIDs();
	output.v2_intime_meanz_v          = selectionv2algo.getIntimeMeanZ();
	output.v2_intime_zfwhm_v          = selectionv2algo.getIntimeZFWHM();
	output.v2_intime_pemax_v          = selectionv2algo.getIntimePEMax();
	output.v2_track_zdiff_frac_v      = selectionv2algo.getIntimeZDiffFrac();

	output.num_matchable_flashes      = selectionv2algo.getMatchableFlashes();
	output.num_matched_flashes        = selectionv2algo.getMatchedFlashes();
      }
    }
    
    
    // ------------------------------------------------------------------------//
    // Make Combined Tagged Image

    if ( m_config.croi_write_cfg.get<bool>("WriteCombinedTaggedImage") ) {

      // create blank image for combined tags
      for ( size_t p=0; p<input.img_v.size(); p++ ) {
        larcv::Image2D combined( input.img_v.at(p).meta() );
        combined.paint(0.0);
        output.combined_v.emplace_back( std::move(combined) );
      }

      for ( size_t itrack=0; itrack<output.flashdata_v.size(); itrack++ ) {
        const larlitecv::TaggerFlashMatchData& flashdata = output.flashdata_v.at(itrack);
        if ( flashdata.m_track3d.NumberTrajectoryPoints()==0 ) {
          continue;
        }
        const std::vector<larcv::Pixel2DCluster>& pixels = flashdata.m_pixels;
        int tagval = 10.0*((int)flashdata.m_type + 1); // 0= nothing, 10 = thrumu, 20=stopmu, 30=untagged/contained, 40=selected CROI
        if ( output.flashdata_selected_v.at(itrack)==1 ) {
          tagval = 40.0;
        }
        for ( size_t p=0; p<output.combined_v.size(); p++ ) {
          const larcv::ImageMeta& meta = input.img_v.at(p).meta();
          for ( auto const& pix : pixels.at(p) ) {
            if ( (int)pix.X()<0 ||  pix.X()>=meta.cols() || (int)pix.Y()<0 || pix.Y()>=meta.rows() ) continue;
            int pixval = output.combined_v.at(p).pixel( pix.Y(), pix.X() );
            if ( tagval > pixval )
              output.combined_v.at(p).set_pixel( pix.Y(), pix.X(), tagval );
          }
        }
      }
    }

    if ( m_config.croi_write_cfg.get<bool>("WriteTrackOpFlashes")) {
      // collect flash info
      output.track_opflash_v.clear();
      for ( auto const& ophypo : selectionalgo.getOpFlashHypotheses() )
        output.track_opflash_v.push_back( ophypo );
    }

    if ( m_config.verbosity>0 )
      std::cout << "== End of Untagged Cluster and CROI Selection ===============================" << std::endl;

    return output;
  }

  void TaggerCROIAlgo::printTimeTracker( int num_events ) {
    const std::string stage_names[] = { "ThruMuConfig", "ThruMuContour", "ThruMuBMT", "ThruMuFlash", "ThruMuFilter", "ThruMuTracker", "StopMuTracker", "Recluster", "Untagged", "PCAmerge", "CROI" };
    float tot_time = 0.;
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "TaggerCROIAlgoConfig::printTimeTracker" << std::endl;
    std::cout << "Number of Events: " << num_events << std::endl;
    for (int i=0; i<kNumStages; i++) {
      std::cout << stage_names[i] << " : " << m_time_tracker[i] << " secs";
      if ( num_events>0 )
        std::cout << "  " << m_time_tracker[i]/float(num_events) << " secs/event";
      std::cout << std::endl;
      tot_time += m_time_tracker[i];
    }
    std::cout << "Total: " << tot_time << " secs";
    if ( num_events>0 )
      std::cout << "  " << tot_time/float(num_events) << " secs/event";
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
  }
}
