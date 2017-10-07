#include "CosmicRetagger.h"

// larlite
#include "DataFormat/track.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "DataFormat/Pixel2DCluster.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

// larlitecv
#include "TaggerTypes/Path2Pixels.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "ThruMu/Pixel2SpacePoint.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "UnipolarHack/UnipolarHackAlgo.h"
#include "TaggerCROIAlgo.h"

namespace larlitecv {

  CosmicRetagger::CosmicRetagger( std::string retagger_cfg, std::string taggerout_larcv_file, std::string taggerout_larlite_file )
    : m_last_loaded_entry(-1)
    , m_last_loaded_image_result(false)
  {
    setConfigFile( retagger_cfg );

    // we have a custom IOManager configuration. It reads every tree as we don't know what trees were written by the version of the tagger un.
    //std::string larcv_io_cfg=" Verbosity: 2 IOMode: 0 InputFiles: [] ReadOnlyDataType: [] ReadOnlyDataName: []";
    //std::string larlite_io_cfg=" Verbosity: 2 IOMode: 0 ReadOnlyProducers: [] ReadOnlyDataTypes: []";

    std::string larcv_io_cfg=" Verbosity: 2 IOMode: 0";
    std::string larlite_io_cfg=" Verbosity: 2 IOMode: 0";
    
    larcv::PSet larcv_pset   = larcv::PSet( "IOManager", larcv_io_cfg );
    larcv::PSet larlite_pset = larcv::PSet( "StorageManager", larlite_io_cfg );

    // input
    m_dataco_input.configure( larcv_pset, larlite_pset );
    m_dataco_input.add_inputfile( taggerout_larcv_file, "larcv" );
    m_dataco_input.add_inputfile( taggerout_larlite_file, "larlite" );
    m_dataco_input.initialize();

    // output
    m_dataco_output.configure( m_cfg_file, "StorageManagerOut", "IOManagerOut", "TaggerCROI" );
    m_dataco_output.initialize();

    m_nentries = m_dataco_input.get_nentries("larcv");
    std::cout << "CosmicRetagger Entries: larcv=" << m_nentries << " larlite=" << m_dataco_input.get_nentries("larlite") << std::endl;

  }

  bool CosmicRetagger::processInputImages() {
    readData();
    readInputData();
    readThrumuData();
    return true;
  }

  bool CosmicRetagger::readData() {
    
    // we look into the file and check what's inside
    if ( !m_state.configured ) {
      std::cout << "Invalid state to run processInputImages: " << printState() << std::endl;
      m_last_loaded_image_result = false;      
      return false;
    }
        
    // -----------------------------------------------------------------------------
    // Load the data
    if ( m_last_loaded_entry==m_entry ) {
      std::cout << "Already Loaded Entry " << m_entry << std::endl;
      return m_last_loaded_entry;
    }
    
    int run, subrun, event;
    try {
      m_dataco_input.goto_entry( m_entry, "larcv" );

      // set the id's
      m_dataco_input.get_id(run,subrun,event);
      m_dataco_output.set_id(run,subrun,event);
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading the input data: " << e.what() << std::endl;
      m_last_loaded_image_result = false;
      return false;
    }
      
    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << m_entry << " : " << run << " " << subrun << " " << event << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // set the state: we reset all downstream states
    m_state.input_ready = m_state.thrumu_run = m_state.stopmu_run = m_state.untagged_run = m_state.croi_run = false;
    return true;
  }
  
  bool CosmicRetagger::readInputData() {

    if ( m_entry<0 ) {
      std::cout << "Entry not loaded from disk yet." << std::endl;
      return false;
    }
    
    // LOAD INPUT_DATA
    m_input_data.clear();

    bool hastpc   = false;
    bool hasbadch = false;
    
    // need a meta to reinterpret some info
    larcv::ImageMeta blank;
    m_meta = blank;

    // get images (from larcv)    
    larcv::EventImage2D* event_imgs   = NULL;
    try {
      event_imgs = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "modimg" );
    }
    catch (...) {
      std::cout << "could not load modimg" << std::endl;
      event_imgs = NULL;
    }
      
    if ( event_imgs!=NULL ) {
      event_imgs->Move( m_input_data.img_v );
      hastpc = true;
    }
    if ( m_input_data.img_v.size()==3 )
      m_meta = m_input_data.img_v.front().meta();

    // badch/gapch channels
    larcv::EventImage2D* event_gapchs = NULL;
    try {
      event_gapchs = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "gapchs" );
    }
    catch (...) {
      std::cout << "could not load gapchs" << std::endl;
      event_gapchs = NULL;
    }
    if ( event_gapchs!=NULL ) {
      event_gapchs->Move( m_input_data.gapch_v );
      for ( auto const& img : m_input_data.gapch_v ) {
	m_input_data.badch_v.push_back( img );
      }
      hasbadch = true;
    } else if ( m_meta.rows()!=0 ) {
      // make empty channels
      for (size_t p=0; p>m_input_data.img_v.size(); p++) {
	larcv::Image2D dummy( m_input_data.img_v[p].meta() );
	m_input_data.gapch_v.push_back( dummy );
	m_input_data.badch_v.push_back( dummy );
      }
    }
    if ( m_meta.rows()==0 && m_input_data.gapch_v.size()==3 ) {
      // if meta not define yet, but we got images, define the meta
      m_meta = m_input_data.gapch_v.front().meta();
    }

    if ( hastpc && hasbadch ) {
      m_state.input_ready = true;
    }

    // opflash
    larlite::event_opflash* pevent_opflash_v = (larlite::event_opflash*)m_dataco_input.get_larlite_data( larlite::data::kOpFlash, "opflash" );
    if ( pevent_opflash_v ) 
      m_input_data.opflashes_v.push_back( pevent_opflash_v );

    return true;
  }

  bool CosmicRetagger::readThrumuData() {
    
   // reload ThruMu payload data: need track points and pixels. and tagged image.
    m_thrumu_data.clear();

    bool has_thrumu_unusedspts  = false;
    bool has_thrumu_tagged      = false;

    // tagged images for boundary point finding (usually not saved)
    larcv::EventImage2D* ev_boundarypixel_v = NULL;
    try {
      ev_boundarypixel_v = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "boundarypixels" );
    }
    catch (...) {
      ev_boundarypixel_v = NULL;
    }
    if ( ev_boundarypixel_v ) {
      ev_boundarypixel_v->Move( m_thrumu_data.boundarypixel_image_v );
    }
    
    larcv::EventImage2D* ev_realspacehits_v = NULL;
    try {
      ev_realspacehits_v = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "realspacehits" );
    }
    catch (...) {
      ev_realspacehits_v = NULL;
    }
    if ( ev_realspacehits_v ) {
      ev_realspacehits_v->Move( m_thrumu_data.realspacehit_image_v );
    }
    
    // load the saved space point

    // prefiltered space points
    larcv::EventPixel2D* prefilter_endpts = NULL;
    try {
      prefilter_endpts = (larcv::EventPixel2D*)m_dataco_input.get_larcv_data( larcv::kProductPixel2D, "prefilterpts" );
    }
    catch (...) {
      prefilter_endpts = NULL;
    }
    if ( prefilter_endpts ) {
      size_t pixplanes = prefilter_endpts->Pixel2DArray().size();
      size_t npts = prefilter_endpts->Pixel2DArray(0).size();
      for (size_t ipt=0; ipt<npts; ipt++) {
	std::vector<larcv::Pixel2D> pixels;
	for (size_t p=0; p<pixplanes; p++) {
	  pixels.push_back(prefilter_endpts->Pixel2DArray(p).at(ipt)); 
	}
	larlitecv::BoundaryEnd_t endpt_t = (larlitecv::BoundaryEnd_t) int(pixels.front().Intensity());
	larlitecv::BoundarySpacePoint sp = larlitecv::Pixel2SpacePoint( pixels, endpt_t, m_input_data.img_v.front().meta() );
	if ( endpt_t<=larlitecv::kDownstream )
	  m_thrumu_data.side_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kAnode )
	  m_thrumu_data.anode_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kCathode )
	  m_thrumu_data.cathode_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kImageEnd )
	  m_thrumu_data.imgends_spacepoint_v.push_back( sp );
      }//end of pt loop
    }//if prefilter container exists

    
    larcv::EventPixel2D* postfilter_endpts[larlitecv::kNumEndTypes];
    std::string endpt_prod_names[7] = { "topspacepts",
					"botspacepts",
					"upspacepts",
					"downspacepts",
					"anodepts",
					"cathodepts",
					"imgendpts" };
    
    for (int endpt_t=0; endpt_t<(int)larlitecv::kNumEndTypes; endpt_t++) {
      postfilter_endpts[endpt_t] = NULL;
      try {
	postfilter_endpts[endpt_t] = (larcv::EventPixel2D*)m_dataco_input.get_larcv_data( larcv::kProductPixel2D, endpt_prod_names[endpt_t] );
      }
      catch (...) {
	postfilter_endpts[endpt_t] = NULL;
      }
      if ( postfilter_endpts[endpt_t]==NULL )
	continue;
      int nendpts = postfilter_endpts[endpt_t]->Pixel2DArray(0).size();
      for (int ipt=0; ipt<nendpts; ipt++) {
	std::vector< larcv::Pixel2D > pixels;
	for (int p=0; p<(int)postfilter_endpts[endpt_t]->Pixel2DArray().size(); p++) {
	  pixels.push_back( postfilter_endpts[endpt_t]->Pixel2DArray(p).at(ipt) );
	}
	if ( m_meta.rows()==0 ) {
	  // set this if the information wasn't available before
	  try {
	    m_meta = postfilter_endpts[endpt_t]->Meta(0);
	  }
	  catch (...){
	  }
	}
	
	std::cout << endpt_prod_names[endpt_t] << " #" << ipt << ": "
		  << " (" << pixels[0].X() << "," << pixels[0].Y() << ")"
		  << " (" << pixels[1].X() << "," << pixels[1].Y() << ")"
		  << " (" << pixels[2].X() << "," << pixels[2].Y() << ")"
		  << std::endl;
	  
	larlitecv::BoundarySpacePoint sp = larlitecv::Pixel2SpacePoint( pixels, (larlitecv::BoundaryEnd_t)endpt_t, m_meta );
	if ( (larlitecv::BoundaryEnd_t)endpt_t<=larlitecv::kDownstream )
	  m_thrumu_data.side_filtered_v.push_back( sp );
	else if ( endpt_t==larlitecv::kAnode )
	  m_thrumu_data.anode_filtered_v.push_back( sp );
	else if ( endpt_t==larlitecv::kCathode )
	  m_thrumu_data.cathode_filtered_v.push_back( sp );
	else if ( endpt_t==larlitecv::kImageEnd )
	  m_thrumu_data.imgends_filtered_v.push_back( sp );
	else
	  continue;
      }
    }

    // Fill EndPts
    larcv::EventPixel2D* unused_endpts[larlitecv::kNumEndTypes];
    for (int endpt_t=0; endpt_t<(int)larlitecv::kNumEndTypes; endpt_t++) {
      unused_endpts[endpt_t] = NULL;
      try {
	unused_endpts[endpt_t] = (larcv::EventPixel2D*)m_dataco_input.get_larcv_data( larcv::kProductPixel2D, std::string("unused_"+endpt_prod_names[endpt_t]) );
      }
      catch (...) {
	unused_endpts[endpt_t] = NULL;
      }
      if ( unused_endpts[endpt_t]==NULL )
	continue;      
      int nendpts = unused_endpts[endpt_t]->Pixel2DArray(0).size();
      for (int ipt=0; ipt<nendpts; ipt++) {
	std::vector< larcv::Pixel2D > pixels;
	for (int p=0; p<(int)unused_endpts[endpt_t]->Pixel2DArray().size(); p++) {
	  pixels.push_back( unused_endpts[endpt_t]->Pixel2DArray(p).at(ipt) );
	}
	if ( m_meta.rows()==0 ) {
	  // set this if the information wasn't available before
	  try {
	    m_meta = unused_endpts[endpt_t]->Meta(0);
	  }
	  catch (...){
	  }	  
	}
	larlitecv::BoundarySpacePoint sp = larlitecv::Pixel2SpacePoint( pixels, (larlitecv::BoundaryEnd_t)endpt_t, m_meta );	
	if ( (larlitecv::BoundaryEnd_t)endpt_t<=larlitecv::kDownstream )
	  m_thrumu_data.side_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kAnode )
	  m_thrumu_data.anode_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kCathode )
	  m_thrumu_data.cathode_spacepoint_v.push_back( sp );
	else if ( endpt_t==larlitecv::kImageEnd )
	  m_thrumu_data.imgends_spacepoint_v.push_back( sp );
	else
	  continue;	
	m_thrumu_data.unused_spacepoint_v.push_back( sp );
      }
      if ( m_thrumu_data.unused_spacepoint_v.size()>0 )
	has_thrumu_unusedspts = true;
    }
    
    larlite::event_track* ev_thrumu_tracks = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    if ( ev_thrumu_tracks!=NULL ) {
      for (size_t itrack=0; itrack<ev_thrumu_tracks->size(); itrack++) {
	m_thrumu_data.track_v.push_back( ev_thrumu_tracks->at(itrack) );
      }
    }

    larcv::EventPixel2D*  ev_thrumu_pixels = NULL;
    try {
      ev_thrumu_pixels = (larcv::EventPixel2D*)m_dataco_input.get_larcv_data( larcv::kProductPixel2D, "thrumupixels" );
    }
    catch (...) {
      ev_thrumu_pixels = NULL;
    }
    if ( ev_thrumu_pixels!=NULL ) {
      for (size_t itrack=0; itrack<ev_thrumu_pixels->Pixel2DClusterArray(0).size(); itrack++) {
	std::vector< larcv::Pixel2DCluster > pixel_v;
	for (int p=0; p<3; p++) {
	  pixel_v.push_back( ev_thrumu_pixels->Pixel2DClusterArray(p).at(itrack) );
	}
	m_thrumu_data.pixelcluster_v.emplace_back( pixel_v );
      }
    }

    larcv::EventImage2D* ev_thrumu_marked = NULL;
    try {
      ev_thrumu_marked = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "marked3d" );
    }
    catch (...) {
      ev_thrumu_marked = NULL;
    }
    if ( ev_thrumu_marked!=NULL ) {
      ev_thrumu_marked->Move( m_thrumu_data.tagged_v );
      has_thrumu_tagged = true;      
    }//end of plane loop
    else if ( m_meta.rows()!=0 ) {
      // we have a meta to make images
      for (int p=0; p<3; p++) {
	larcv::Image2D marked( m_meta );
	marked.paint(0);
	m_thrumu_data.tagged_v.emplace_back( marked );
      }
      if ( m_thrumu_data.pixelcluster_v.size()>0 ) {      
	// can make a tagged image using saved pixel clusters      	
	for (int itrack=0; itrack<(int)m_thrumu_data.pixelcluster_v.size(); itrack++) {
	  for (int p=0; p<3; p++) {
	    const larcv::Pixel2DCluster& pixels = m_thrumu_data.pixelcluster_v.at(itrack).at(p);
	    for ( auto const& pix : pixels ) {
	      m_thrumu_data.tagged_v[p].set_pixel( pix.Y(), pix.X(), 10.0 );
	    }
	  }
	}
	has_thrumu_tagged = true;
      }
      else if ( m_thrumu_data.track_v.size()>0 && m_input_data.img_v.size()>0) {
	// make tagged image from larlite tracks
	bool also_fill_pixel_v = false;
	if ( m_thrumu_data.pixelcluster_v.size()==0 )
	  also_fill_pixel_v = true;
	for (int itrack=0; itrack<(int)m_thrumu_data.track_v.size(); itrack++) {
	  const larlite::track& lltrack = m_thrumu_data.track_v[itrack];
	  std::vector< std::vector<double> > path3d( lltrack.NumberTrajectoryPoints());
	  for (size_t n=0; n<lltrack.NumberTrajectoryPoints(); n++) {
	    path3d[n].resize(3,0.0);
	    const TVector3& pos = lltrack.LocationAtPoint(n);
	    for (int i=0; i<3; i++)
	      path3d[n][i] = pos[i];
	  }
	  std::vector<larcv::Pixel2DCluster> pixels = larlitecv::getTrackPixelsFromImages( path3d, m_input_data.img_v, m_input_data.gapch_v,
											   m_tagger_cfg.thrumu_tracker_cfg.pixel_threshold, m_tagger_cfg.thrumu_tracker_cfg.tag_neighborhood,
											   0.3 );
	  if ( also_fill_pixel_v ) {
	    m_thrumu_data.pixelcluster_v.emplace_back( std::move(pixels) );
	  }
	  for (size_t p=0; p<pixels.size(); p++) {
	    for ( auto const& pix : pixels[p] ) {
	      m_thrumu_data.tagged_v[p].set_pixel( pix.Y(), pix.X(), 10.0 );
	    }
	  }
	  has_thrumu_tagged = true;	  
	}//end of track look
      }
    }//end of if meta available to define images

    if ( has_thrumu_unusedspts && has_thrumu_tagged ) {
      // we have the thrumu info needed to run stopmu
      m_state.thrumu_run = true;
    }
    
    return true;    
  }

  bool CosmicRetagger::readCROIData() {
    // CROI
    larlite::event_track* evout_tracks_alltracks      = NULL;
    try {
      evout_tracks_alltracks = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "all3d" );
    }
    catch (...){
    }
    
    larlite::event_track* evout_tracks_mergedstopmu3d = NULL;
    try {
      evout_tracks_mergedstopmu3d = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "mergedstopmu3d" );
    }
    catch (...){
    }
    
    larlite::event_track* evout_tracks_mergedthrumu3d = NULL;
    try {
      evout_tracks_mergedthrumu3d = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "mergedthrumu3d" );
    }
    catch (...){
    }
    
    larlite::event_track* evout_tracks_mergeduntagged = NULL;
    try {
      evout_tracks_mergeduntagged = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "mergeduntagged3d" );
    }
    catch (...){
    }
    
  }
}
