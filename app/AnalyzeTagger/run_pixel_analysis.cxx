/*
  TAGGER PIXEL ANALYSIS

  Inputs:

  larcv source file
  -----------------
  original images. (set in cfg with InputLArCVImages)
  segmentation imgs. (only if IsMC: true )
  chstatus object.
  badch images. (can specify to get from file or make from chstatus images)
  
  larlite source file
  -------------------
  mctruth from "generator"
  mctrack from "mcreco"
  mcshower from "mcreco"
  trigger from "triggersim" or "trigger" (set in cfg with TriggerProducerName)

  larcv tagger/CROI file
  ----------------------
  pixel2dclusters for thrumu tracks "thrumupixels"
  pixel2dclusters for stopmu tracks "stopmupixels"
  pixel2dclusters for untagged clusters "untaggedpixels"
  pixel2dclusters for contained ROI "croipixels"
  ROIs for contained ROIs "croi"
  boundary space points (only relevant if MC info available, IsMC:true)

  larlite tagger/CROI file
  -------------------------
  not used

 */



#include <iostream>
#include <string>
#include <sstream>

// ROOT
#include "TFile.h"
#include "TTree.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventChStatus.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "UBWireTool/UBWireTool.h"

// larlite
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"
#include "DataFormat/trigger.h"
#include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "TaggerTypes/dwall.h"
#include "MCTruthTools/crossingPointsAnaMethods.h"
#include "MCTruthTools/extractTruthMethods.h"



// OpenCV
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif




int main( int nargs, char** argv ) {

  // PIXEL ANALYSIS

  if ( nargs!=2 ) {
    std::cout << "usage: ./run_pixel_analysis_old_version [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("AnalyzeTagger");

  enum SourceTypes_t { kSource=0, kCROIfile, kNumSourceTypes };
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "LOADING " << source_param[isrc] << " FILES" << std::endl;
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "AnalyzeTagger" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }
  
  // =====================================================================
  // configuration parameters

  std::string outfname   = pset.get<std::string>("OutputAnaFile");
  float fthreshold       = pset.get<float>("PixelThreshold");
  int tag_neighborhood   = pset.get<int>("TagNeighborhood",10);
  int fvertex_radius     = pset.get<int>("PixelRadius");
  int verbosity          = pset.get<int>("Verbosity",0);
  bool ismc              = pset.get<bool>("IsMC");
  bool badch_in_file     = pset.get<bool>("BadChImageInFile");
  std::string inputimgs  = pset.get<std::string>("InputLArCVImages");
  std::string trigname   = pset.get<std::string>("TriggerProducerName");
  bool printFlashEnds    = pset.get<bool>("PrintFlashEnds");
  bool use_reclustered   = pset.get<bool>("UseReclustered");  
  std::vector<std::string> flashprod  = pset.get<std::vector<std::string> >("OpFlashProducer");

 // =====================================================================

  // setup input
  int kThruMu, kStopMu, kUntagged, kCROI, kNumStages;

  std::vector<std::string> stages_pixel_producers;
  std::vector<std::string> stages_track_producers;  
  if ( use_reclustered ) {
    kNumStages = 4;
    kThruMu = 0;
    kStopMu = 1;    
    kUntagged = 2;
    kCROI = 3;
    stages_pixel_producers.resize(kNumStages);
    stages_pixel_producers[0] = "mergedthrumupixels";
    stages_pixel_producers[1] = "mergedstopmupixels";
    stages_pixel_producers[2] = "mergeduntaggedpixels";
    stages_pixel_producers[3] = "croipixels";
    stages_track_producers.resize(kNumStages);    
    stages_track_producers[0] = "mergedthrumu3d";
    stages_track_producers[1] = "mergedstopmu3d";
    stages_track_producers[2] = "mergeduntagged3d";
    stages_track_producers[3] = "croi3d";    
  }
  
  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* tree = new TTree("pixana", "Pixel-level analysis");

  // Event Index
  int run, subrun, event;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;

  // Truth Pixel Quantities
  int nthreshold_pixels[4]; // number of pixels above threshold
  int ncosmic_pixels[4];    // number of non-neutrino pixels
  int nnu_pixels[4];        // number of neutrino pixels
  int nvertex_pixels[4];    // number of neutrino pixels within some pixel radius of vertex
  int nvertex_badch[4];

  // Tagged Pixel Quantities
  int ncosmic_tagged[kNumStages][4];       // number of non-neutrino pixels tagged
  int ncosmic_tagged_once[kNumStages][4];  // number of non-neutrino pixels tagged
  int ncosmic_tagged_many[kNumStages][4];  // number of non-neutrino pixels tagged
  int nnu_tagged[kNumStages][4];           // number of neutrino pixels tagged
  int nvertex_tagged[kNumStages][4];       // number of neutrino pixels within some pixel radius of vertex tagged
  int nvertex_incroi[4];                   // number of pixels near neutrino vertex that are in an ROI
  std::stringstream s_arr;
  s_arr << "[" << (int)kNumStages << "][4]/I";

  // ROI quantities
  int num_rois;     // number of identified ROis
  int nnu_inroi[4]; // number of nu pixels contained in the CROI
  int vertex_in_croi; // is vertex in an CROI
  float closest_dist_to_vertex;
  int closest_dist_stage;

  // Crossing Point data
  larlitecv::CrossingPointAnaData_t xingptdata;

  // Event
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");

  // Truth
  truthdata.bindToTree( tree );

  tree->Branch("nthreshold_pixels", nthreshold_pixels, "nthreshold_pixels[4]/I");
  tree->Branch("ncosmic_pixels",    ncosmic_pixels,    "ncosmic_pixels[4]/I");
  tree->Branch("nnu_pixels",        nnu_pixels,        "nnu_pixels[4]/I");
  tree->Branch("nvertex_pixels",    nvertex_pixels,    "nvertex_pixels[4]/I");
  tree->Branch("nvertex_badch",     nvertex_badch,     "nvertex_badch[4]/I");

  tree->Branch("ncosmic_tagged",      ncosmic_tagged,      std::string("ncosmic_tagged"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_once", ncosmic_tagged_once, std::string("ncosmic_tagged_once"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_many", ncosmic_tagged_many, std::string("ncosmic_tagged_many"+s_arr.str()).c_str() );
  tree->Branch("nnu_tagged",          nnu_tagged,          std::string("nnu_tagged"+s_arr.str()).c_str() );
  tree->Branch("nvertex_tagged",      nvertex_tagged,      std::string("nvertex_tagged"+s_arr.str()).c_str() );
  tree->Branch("nvertex_incroi",      nvertex_incroi,      "nvertex_incroi[4]/I" );

  tree->Branch("num_rois", &num_rois, "num_rois/I");
  tree->Branch("nnu_inroi", nnu_inroi, "nnu_inroi[4]" );
  tree->Branch("vtx_in_croi", &vertex_in_croi,         "vtx_in_croi/I" );
  tree->Branch("dist_to_vtx", &closest_dist_to_vertex, "dist_to_vtx/F" );
  tree->Branch("stage_at_vtx", &closest_dist_stage, "stage_at_vtx/I" );
  
  xingptdata.bindToTree( tree );
  
  // Space Charge Corrections
  larlitecv::SpaceChargeMicroBooNE sce;

  // Empty Channel Algo
  larlitecv::EmptyChannelAlgo emptyalgo;

  int nentries = dataco[kCROIfile].get_nentries("larcv");

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
  if ( startentry>=nentries ) {
    std::cout << "Starting beyond end of file. Nothing to do." << std::endl;
    return 0;
  }
  std::cout << "Start Entry: " << startentry << std::endl;
  std::cout << "End Entry: " << endentry-1 << std::endl;
  std::cout << "Buckle up!" << std::endl;
  

  for (int ientry=startentry; ientry<endentry; ientry++) {

    // load the entry
    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(run,subrun,event);

    // sync up the other files
    for ( int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
      if ( isrc!=kCROIfile ) {
        dataco[isrc].goto_event(run,subrun,event, "larcv");
      }
    }

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << ": (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    // initialize the output variables
    truthdata.clear();
    xingptdata.clear();
    for (int p=0; p<4; p++) {
      ncosmic_pixels[p] = 0;
      nnu_pixels[0] = 0;
      nvertex_pixels[p] = 0;
      nthreshold_pixels[p] = 0;
      nnu_inroi[p] = 0;
      nvertex_badch[p] = 0;
      for (int istage=0; istage<kNumStages; istage++) {
	ncosmic_tagged[istage][p] = 0;
	ncosmic_tagged_once[istage][p] = 0;
	ncosmic_tagged_many[istage][p] = 0;
	nnu_tagged[istage][p] = 0;
	nvertex_tagged[istage][p] = 0;
      }
      nvertex_incroi[p] = 0;
    }
    vertex_in_croi = 0;    


    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,inputimgs);
    larcv::EventImage2D* ev_badch  = NULL;
    larcv::EventImage2D* ev_segs   = NULL;
    if ( ismc ) {
      ev_segs = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    }
    if ( badch_in_file ) {
      // we either get the badch from the file
      ev_badch = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"gapchs");
    }
    else {
      // or we have to make the badch image from a ChStatus object
      larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "tpc" );
      std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      ev_badch = new larcv::EventImage2D;
      ev_badch->Emplace( std::move(chstatus_img_v) );
    }

    // get the output of the tagger
    larcv::EventPixel2D* ev_pix[kNumStages] = {0};
    larlite::event_track* ev_track[kNumStages] = {0};
    try {
      ev_pix[kThruMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kThruMu]);
      ev_track[kThruMu]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kThruMu]);      
    }
    catch (...) {
      ev_pix[kThruMu] = NULL;
      ev_track[kThruMu] = NULL;
    }
    try {
      ev_pix[kStopMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kStopMu]);
      ev_track[kStopMu]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kStopMu]);            
    }
    catch (...) {
      ev_pix[kStopMu] = NULL;
    }
    try {
      ev_pix[kUntagged]  = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kUntagged]);
      ev_track[kUntagged]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kUntagged]);                  
    }
    catch (...) {
      ev_pix[kUntagged] = NULL;
    }
    try {
      ev_pix[kCROI]      = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kCROI]);
      ev_track[kCROI]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kCROI]);                  
    }
    catch (...) {
      ev_pix[kCROI] = NULL;
    }

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& chstatus_img_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>* segs_v   = NULL;
    if ( ismc ) {
      segs_v = &(ev_segs->Image2DArray());
    }

    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = NULL;
    try {
      ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    }
    catch (...) {
      ev_contained_roi = NULL;
    }
    std::vector<larcv::ROI> containedrois_v;
    if ( ev_contained_roi ) {
      containedrois_v = ev_contained_roi->ROIArray();
      num_rois = (int)containedrois_v.size();
      std::cout << "====ROIs===========================" << std::endl;
      for ( auto& roi : containedrois_v ) {
	std::cout << " roi: " << roi.dump();
      }
      std::cout << "===================================" << std::endl;
    }
      
    // get the boundary end point info (only if have MC info to compare against)
    std::vector<larcv::EventPixel2D*> ev_spacepoints(7,0);
    std::string spacepoint_producers[7] = { "topspacepts", "botspacepts", "upspacepts", "downspacepts", "anodepts", "cathodepts", "imgendpts" };
    for ( int i=0; i<7; i++ ) {
      try {
	ev_spacepoints[i] = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,spacepoint_producers[i]);
      }
      catch (...) {
	ev_spacepoints[i] = NULL;
      }
      if ( ev_spacepoints[i]!=NULL )
	std::cout << "number of " << spacepoint_producers[i] << ": " << ev_spacepoints[i]->Pixel2DArray(0).size() << std::endl;	
    }

    // get the opflashes
    std::vector< larlite::event_opflash* > opflash_v;
    for ( auto const& prodname : flashprod ) {
      larlite::event_opflash* ev_flash = (larlite::event_opflash*)dataco[kSource].get_larlite_data(larlite::data::kOpFlash, prodname);
      std::cout << "number of flashes in " << prodname << ": " << ev_flash->size() << std::endl;
      opflash_v.push_back( ev_flash );
    }

    // get other information, e.g. truth
    larlite::event_mctruth* ev_mctruth   = NULL;
    larlite::event_mctrack* ev_mctrack   = NULL;
    larlite::event_mcshower* ev_mcshower = NULL;
    larlite::trigger* ev_trigger         = (larlite::trigger*) dataco[kSource].get_larlite_data(larlite::data::kTrigger,trigname);
    if ( ismc ) {
      ev_mctruth = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
      ev_mctrack = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
      ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");

      // extract the truth quantities of interest
      larlitecv::extractTruth( *ev_mctruth, *ev_mctrack, truthdata );
    }

    // =======================================================================
    // MC INFO ANALYSIS

    // quantities filled if MC present
    std::vector<larcv::Image2D> nupix_imgs_v;
    // std::vector< std::vector<int> > start_pixels;
    // std::vector< std::vector<float> > start_crossingpts;
    // std::vector< std::vector<int> > end_pixels;
    // std::vector< std::vector<float> > end_crossingpts;
    std::vector<int> vertex_col(3,-1);
    std::vector<double> vtx_sce(3,0);    
    int vertex_row = -1;
      
    if ( ismc ) {
    
      // get the vertex in the pixel coordinates
      std::vector<int> wid(3,-1);
      std::vector<double> dpos(3);
      for (int i=0; i<3; i++ ) dpos[i] = truthdata.pos[i];
      std::vector<double> vtx_offset = sce.GetPosOffsets( dpos[0], dpos[1], dpos[2] );

      for (int i=1; i<3; i++ ) vtx_sce[i] = dpos[i] + vtx_offset[i];
      vtx_sce[0] = dpos[0] - vtx_offset[0] + 0.7;
      
      for (size_t p=0; p<3; p++) {
	wid[p] = ::larutil::Geometry::GetME()->WireCoordinate( vtx_sce, p );
	if ( wid[p]>=0 && wid[p]<3456 )
	  vertex_col[p] = imgs_v.at(p).meta().col(wid[p]);
      }
      float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
      float vertex_tick = vtx_sce[0]/cm_per_tick + 3200.0;
      if ( vertex_tick >= imgs_v.at(0).meta().min_y() && vertex_tick<=imgs_v.at(0).meta().max_y() )
	vertex_row = imgs_v.at(0).meta().row( vertex_tick );

      std::cout << "Vertex Pixel Coordinates (SCE corrected): (" << vertex_row << ", " << vertex_col[0] << "," << vertex_col[1] << "," << vertex_col[2] << ")" << std::endl;

      // did any of the ROIs contain the vertex?
      vertex_in_croi = 0;
      for ( auto& roi : containedrois_v ) {
	int nplanes_in_roi = 0;
	for (size_t p=0; p<3; p++ ) {
	  const larcv::ImageMeta& bb = roi.BB( (larcv::PlaneID_t)p );
	  if ( vertex_tick>=bb.min_y() && vertex_tick<=bb.max_y() && vertex_col[p]>=bb.min_x() && vertex_col[p]<=bb.max_x() )
	    nplanes_in_roi++;
	}
	if (nplanes_in_roi>=2) {
	  vertex_in_croi = 1;
	  break;
	}
      }

      // loop over MC tracks, get end points of muons
      larlitecv::analyzeCrossingMCTracks( xingptdata, imgs_v.front().meta(), imgs_v, ev_trigger, ev_mctrack, opflash_v, printFlashEnds );
      // int intime_cosmics = xingptdata.true_intime_thrumu + xingptdata.true_intime_stopmu;
      // std::cout << "number of intime cosmics: "       << intime_cosmics << std::endl;
      // std::cout << "number of intime thrumu: "        << xingptdata.true_intime_thrumu << std::endl;
      // std::cout << "number of intime stopmu: "        << xingptdata.true_intime_stopmu << std::endl;
      std::cout << "number of true crossing points: " << xingptdata.tot_true_crossingpoints << std::endl;
      for (int i=0; i<6; i++) {
	std::cout << "  " << spacepoint_producers[i] << ": " << xingptdata.true_crossingpoints[i] << std::endl;
      }

      // make truth pixel counts
      // count the pixels. determine if cosmic and neutrino are tagged. also if neutrino is in rois
      // we loop through the rows and cols
      for (size_t p=0; p<3; p++) {
	  
	// we create a neutrino pixel image, to make things easier downstream
	larcv::Image2D nupix_img( imgs_v.at(p).meta() );
	nupix_img.paint(0);
	for (size_t row=0; row<imgs_v.at(p).meta().rows(); row++) {
	  for (size_t col=0; col<imgs_v.at(p).meta().cols(); col++) {
	      
	    // check if this is a pixel of interest
	    if ( imgs_v.at(p).pixel(row,col)<fthreshold ) continue;
	      
	    bool near_vertex = false;
	    bool is_nu_pix = false;
	      
	    // are we some radius from the vertex?
	    if ( (int)row>=vertex_row-fvertex_radius && (int)row<=vertex_row+fvertex_radius
		 && (int)col>=vertex_col[p]-fvertex_radius && (int)col<=vertex_col[p]+fvertex_radius ) {
	      near_vertex = true;
	      nvertex_pixels[p]++;
	      nvertex_pixels[3]++;
	    }
	      
	    // above threshold. is it a neutrino pixel?
	    bool in_seg_image = false;
	    int seg_row = -1;
	    int seg_col = -1;
	      
	    const larcv::Image2D& segimg = segs_v->at(p);
	    float x = imgs_v.at(p).meta().pos_x(col);
	    float y = imgs_v.at(p).meta().pos_y(row);
	    if ( x>segs_v->at(p).meta().min_x() && x<segs_v->at(p).meta().max_x()
		 && y>segs_v->at(p).meta().min_y() && y<segs_v->at(p).meta().max_y() ) {
	      in_seg_image = true;
	      seg_row = segs_v->at(p).meta().row(y);
	      seg_col = segs_v->at(p).meta().col(x);
	    }

	    if ( in_seg_image && segs_v->at(p).pixel(seg_row,seg_col)>0 ) {
		
	      is_nu_pix = true;
	      nnu_pixels[p]++;
	      nnu_pixels[3]++;

	      // is the neutrino pixel inside the ROI?
	      for ( auto const& cand_roi : containedrois_v ) {
		float wired = imgs_v.at(p).meta().pos_x(col);
		float tick  = imgs_v.at(p).meta().pos_y(row);
		const larcv::ImageMeta& cand_roi_bb = cand_roi.BB().at(p);
		if ( cand_roi_bb.min_x()<wired && wired<cand_roi_bb.max_x()
		     && cand_roi_bb.min_y()<tick && tick<cand_roi_bb.max_y() )
		  {
		    nnu_inroi[p]++;
		    nnu_inroi[3]++;
		  }
	      }
	    }//if in semgment image
	    else {
	      // not a neutrino, so cosmic
	      ncosmic_pixels[p]++;
	      ncosmic_pixels[3]++;
	    }//end if cosmic
	      
	    if ( is_nu_pix )
	      nupix_img.set_pixel( row, col, 1.0 );
	    if ( near_vertex )
	      nupix_img.set_pixel( row, col, 10.0 );
	      
	  }//end of col loop
	}//end of row loop
	nupix_imgs_v.emplace_back( std::move(nupix_img) );
      }//end of loop over planes for counting neutrino/cosmic pixels
    }//end of MC Pixel Counts
      
    // did any of the ROIs catch neutrino pixels?
    for ( int r=vertex_row-fvertex_radius; r<=vertex_row+fvertex_radius; r++ ) {
      if ( r<0 || r>=(int)imgs_v.front().meta().rows() ) continue;
      for (size_t p=0; p<3; p++) {
	const larcv::ImageMeta& meta = imgs_v.at(p).meta();
	for ( int c=vertex_col[p]-fvertex_radius; c<=vertex_col[p]+fvertex_radius; c++ ) {
	  if ( c<0 || c>=(int)meta.cols() ) continue;
	  
	  // is this a bad ch pixel?
	  if ( r==vertex_row && chstatus_img_v.at(p).pixel(r,c)>0 ) {
	    nvertex_badch[p]++;
	    nvertex_badch[3]++;
	  }
	    
	  if ( imgs_v.at(p).pixel(r,c)>fthreshold ) {
	      
	    float wire = meta.pos_x( c );
	    float tick = meta.pos_y( r );
	      
	    // a vertex pixel!
	    // let's loop over ROIs
	    bool inroi = false;
	    for ( auto const& croi : containedrois_v ) {
	      if ( wire>=croi.BB(p).min_x() && wire<=croi.BB(p).max_x()
		   && tick>=croi.BB(p).min_y() && tick<=croi.BB(p).max_y() ) {
		inroi = true;
		break;
	      }
	    }
	    if ( inroi ) {
	      nvertex_incroi[p]++;
	      nvertex_incroi[3]++;
	    }
	  }
	}//end of col loop
      }//end of plane lopp
    }//end of vertex row loop

    std::cout << "Number of vertex pixels in an ROI: " << nvertex_incroi[3] << " out of " << nvertex_pixels[3] << std::endl;
    std::cout << "Fraction of vertex pixels are badchannels: " << float(nvertex_badch[3])/float(nvertex_pixels[3]) << std::endl;

    analyzeCrossingDataOnly( xingptdata, ev_spacepoints );
    analyzeCrossingMatches( xingptdata,  ev_spacepoints, imgs_v.front().meta() );
    
    // ==========================================================================================


    // non-MC Pixel counting
    for (size_t p=0; p<3; p++) {
      for (size_t row=0; row<imgs_v.at(p).meta().rows(); row++) {
	for (size_t col=0; col<imgs_v.at(p).meta().cols(); col++) {
	  
	  // check if this is a pixel of interest
	  if ( imgs_v.at(p).pixel(row,col)<fthreshold ) continue;
	    
	  nthreshold_pixels[p]++;
	  nthreshold_pixels[3]++;
	  if ( !ismc ) {
	    // if no MC, every pixel is a cosmic pixel
	    ncosmic_pixels[p]++;
	    ncosmic_pixels[3]++;	    
	  }
	}
      }
    }//end of plane loop


    // Loop over track pixels
    // make left over image
#ifdef USE_OPENCV
    std::vector<cv::Mat> cvleftover_v;
    for ( auto const& img : imgs_v ) {
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, fthreshold, 100.0 );
      cvleftover_v.emplace_back(std::move(cvimg));
    }
#endif

    // we need an image to mark how times a pixel has been marked
    std::vector< larcv::Image2D > img_marker_v;
    for (size_t p=0; p<3; p++) {
      larcv::Image2D img_counter( imgs_v.at(p).meta() );
      img_counter.paint(0);
      img_marker_v.emplace_back( std::move(img_counter) );
    }

    // loop over tagger stage results
    for ( int istage=0; istage<kNumStages; istage++ ) {
      if ( ev_pix[istage]==NULL )
	continue;
      for ( size_t p=0; p<3; p++ ) {
        // we need images to track the number of times pixels are Tagged
        larcv::Image2D& img_counter = img_marker_v.at(p);
        if ( istage==kCROI ) {
          // we reset the tracker for the CROI pixel counting
          img_counter.paint(0);
        }

        for ( auto const& pixcluster : ev_pix[istage]->Pixel2DClusterArray(p) ) {
	  
          // for each track cluster, we loop over pixels and store unique pixels by using a set
          std::set< std::vector<int> > tagged_set;
          for ( auto const& pix : pixcluster ) {
	    
            for ( int dr=-tag_neighborhood; dr<=tag_neighborhood; dr++ ) {
              int r = pix.Y()+dr;
              if ( r<0 || r>=(int)imgs_v.at(p).meta().rows() ) continue;
              for ( int dc=-tag_neighborhood; dc<=tag_neighborhood; dc++ ) {
                int c = pix.X()+dc;
                if ( c<0 || c>=(int)imgs_v.at(p).meta().cols() ) continue;
                if ( imgs_v.at(p).pixel(r,c)<fthreshold ) continue;
		
                std::vector<int> pixel(2);
                pixel[0] = c;
                pixel[1] = r;
                tagged_set.insert( pixel );
              }
            }
          }//end of pixel loop
          //std::cout << "stage " << istage << " pix in cluster=" << pixcluster.size() << " tagged=" << tagged_set.size() << std::endl;
	  
          // we loop over the set, tagging image.  We increment the counter of each pixel.
          for ( auto const& pix : tagged_set ) {
            float val = img_counter.pixel( pix[1], pix[0]) + 1.0;
            img_counter.set_pixel( pix[1], pix[0], val );
	    
#ifdef USE_OPENCV
            // set color of tagged pixel
            if ( istage==kCROI ) {
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[0] = 255;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[1] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[2] = 255;
            }
            else if ( istage==kUntagged ) {
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[0] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[1] = 255;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[2] = 0;
            }
            else if ( istage==kThruMu ) {
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[0] = 200;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[1] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[2] = 0;
            }
            else if ( istage==kStopMu ) {
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[0] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[1] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[2] = 200;
            }
            else {
              throw std::runtime_error("oops.");
            }
#endif
          }
        }//end of pix cluster loop
	
        // total the pixels
        for ( size_t r=0; r<img_counter.meta().rows(); r++) {
          for ( size_t c=0; c<img_counter.meta().cols(); c++ ) {

	    if ( imgs_v.at(p).pixel(r,c)<fthreshold ) continue;
	    
            // is pixel cosmic or neutrino
	    bool isnupix = false;
	    bool isvertex = false;
	    if ( ismc ) {
	      isnupix  = ( nupix_imgs_v.at(p).pixel(r,c)>=1.0  ) ? true : false;
	      isvertex = ( nupix_imgs_v.at(p).pixel(r,c)>=10.0 ) ? true : false;
	    }
	    
            int count = img_counter.pixel( r, c );
	    
            if ( count>0 ) {
              if ( isnupix ) {
                nnu_tagged[istage][p]++;
                nnu_tagged[istage][3]++;
              }
              if ( isvertex ) {
                nvertex_tagged[istage][p]++;
                nvertex_tagged[istage][3]++;
              }
              if ( !isnupix && !isvertex ) {
                // is cosmic
                ncosmic_tagged[istage][p]++;
                ncosmic_tagged[istage][3]++;
                if (count==1)  {
                  ncosmic_tagged_once[istage][p]++;
                  ncosmic_tagged_once[istage][3]++;
                }
                else if ( count>1 ) {
                  ncosmic_tagged_many[istage][p]++;
                  ncosmic_tagged_many[istage][3]++;
                }
              }
            }
          }//end of loop over c of counter image
        }//loop over rows of counter image
      }//end of plane loop
    }//end of stage loop
    
    // for the totals for the thrumu/stop/untagged, we want the new amount of pixels tagged, not the aggregate past that stage
    for ( int istage=kUntagged; istage>kThruMu; istage-- ) {
      // we save the difference
      for (int p=0; p<4; p++) {
        ncosmic_tagged[istage][p]      = ncosmic_tagged[istage][p] - ncosmic_tagged[istage-1][p];
        ncosmic_tagged_once[istage][p] = ncosmic_tagged_once[istage][p] - ncosmic_tagged_once[istage-1][p];
        ncosmic_tagged_many[istage][p] = ncosmic_tagged_many[istage][p] - ncosmic_tagged_many[istage-1][p];
        nnu_tagged[istage][p]          = nnu_tagged[istage][p] - nnu_tagged[istage-1][p];
        nvertex_tagged[istage][p]      = nvertex_tagged[istage][p] - nvertex_tagged[istage-1][p];
      }
    }


    //================================================================================
    // Analyze if a track came close to a vertex
    if (ismc) {
      closest_dist_to_vertex = 1.0e6;
      closest_dist_stage = -1;
      for (int istage=0; istage<kNumStages; istage++) {
	for (int itrack=0; itrack<(int)(ev_track[istage]->size()); itrack++) {
	  const larlite::track& t = (*ev_track[istage])[itrack];
	  for (int ipt=0; ipt<(int)t.NumberTrajectoryPoints(); ipt++) {
	    float dist=0;
	    for (int v=0; v<3; v++) {
	      float dx = vtx_sce[v]-t.LocationAtPoint(ipt)[v];
	      dist += dx*dx;
	    }
	    dist = sqrt(dist);
	    if ( dist<=closest_dist_to_vertex ) {
	      closest_dist_to_vertex = dist;
	      closest_dist_stage = istage;
	    }
	  }
	}
      }
    }

    //================================================================================

#ifdef USE_OPENCV
    // draw image
    if ( pset.get<bool>("SaveJPEG") ) {
      int color_codes[7][3] = { {255,0,0}, // top
				{0,0,255}, // bot
				{0,255,255}, // upstream
				{0,255,0}, // downstream
				{255,255,0}, // anode
				{255,0,255}, // cathode
				{0,128,255} }; // imageend
      for ( size_t p=0; p<cvleftover_v.size(); p++ ) {
	auto& leftover = cvleftover_v.at(p);

	if ( ismc ) {
	  // draw truth end points!
	  std::cout << "startpixels=" << xingptdata.start_pixels.size() << " " << xingptdata.start_type.size() << std::endl;
	  for ( size_t istart=0; istart<xingptdata.start_type.size(); istart++) {
	    auto const& start_pix = xingptdata.start_pixels[istart];
	    int i = xingptdata.start_type[istart];
	    if ( xingptdata.start_type[istart]<4 )
	      cv::circle( leftover, cv::Point(start_pix[p+1],start_pix[0]), 4, cv::Scalar( color_codes[i][0], color_codes[i][1], color_codes[i][2]), 2, -1 );
	    else
	      cv::rectangle( leftover, cv::Point(start_pix[p+1]-4,start_pix[0]-4), cv::Point(start_pix[p+1]+4,start_pix[0]+4),
			     cv::Scalar( color_codes[i][0], color_codes[i][1], color_codes[i][2]), 1 );
	  }
	  std::cout << "endpixels=" << xingptdata.end_pixels.size() << " " << xingptdata.end_type.size() << std::endl;	
	  for ( size_t iend=0; iend<xingptdata.end_type.size(); iend++) {	
	    auto const& end_pix = xingptdata.end_pixels[iend];
	    int i = xingptdata.end_type[iend];	  
	    if ( xingptdata.end_type[iend]<4 )
	      cv::circle( leftover, cv::Point(end_pix[p+1],end_pix[0]), 4, cv::Scalar( color_codes[i][0], color_codes[i][1], color_codes[i][2]), 2, -1 );
	    else
	      cv::rectangle( leftover, cv::Point(end_pix[p+1]-4,end_pix[0]-4), cv::Point(end_pix[p+1]+4,end_pix[0]+4),
			     cv::Scalar( color_codes[i][0], color_codes[i][1], color_codes[i][2]), 1 );
	  }
	
	  // draw proposed end points
	  for ( int i=0; i<7; i++) {
	    if ( ev_spacepoints[i]==NULL )
	      continue;
	    for ( auto const& endpt : ev_spacepoints[i]->Pixel2DArray(p) ) {
	      cv::drawMarker( leftover, cv::Point(endpt.X(), endpt.Y()),  cv::Scalar( color_codes[i][0], color_codes[i][1], color_codes[i][2]), cv::MARKER_CROSS, 6, 2);
	    }
	  }

	  // sce vertex
	  cv::circle( leftover, cv::Point(vertex_col[p],vertex_row), 4, cv::Scalar(0,0,255),   2, -1 );
	  cv::circle( leftover, cv::Point(vertex_col[p],vertex_row), 3, cv::Scalar(0,255,255), 1, -1 );      

	}

	// draw roi
	for ( auto const& roi : containedrois_v ) {
	  larcv::draw_bb( leftover, imgs_v.front().meta(), roi.BB(p), 255, 0, 255, 2 );
	}

      
	std::stringstream ss;
	ss << "leftover_clust_i" << ientry << "_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
	std::cout << "write: " << ss.str() << std::endl;
	cv::imwrite( ss.str(), leftover );
      }
    }
#endif

    tree->Fill();

  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
