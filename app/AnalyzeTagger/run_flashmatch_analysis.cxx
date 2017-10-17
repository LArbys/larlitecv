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
  boundary space points (only relevant if MC info available, IsMC:true)1

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
#include "Base/DataFormatConstants.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"
#include "DataFormat/trigger.h"
#include "DataFormat/track.h"
#include "DataFormat/user_info.h"
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

  // FLASH MATCH ANALYSIS

  if ( nargs!=2 ) {
    std::cout << "usage: ./run_flashmatch_analysis [config file]" << std::endl;
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

  // setup stage defintions
  // ----------------------
  int kThruMu, kStopMu, kUntagged, kCROI, kAll, kNumStages;

  std::vector<std::string> stages_pixel_producers;
  std::vector<std::string> stages_track_producers;  
  kNumStages = 5;
  kThruMu = 0;
  kStopMu = 1;    
  kUntagged = 2;
  kCROI = 3;
  kAll = 4;  
  stages_pixel_producers.resize(kNumStages);
  stages_pixel_producers[0] = "mergedthrumupixels";
  stages_pixel_producers[1] = "mergedstopmupixels";
  stages_pixel_producers[2] = "mergeduntaggedpixels";
  stages_pixel_producers[3] = "croipixels";
  stages_pixel_producers[4] = "allpixels";
  stages_track_producers.resize(kNumStages);    
  stages_track_producers[0] = "mergedthrumu3d";
  stages_track_producers[1] = "mergedstopmu3d";
  stages_track_producers[2] = "mergeduntagged3d";
  stages_track_producers[3] = "croi3d";
  stages_track_producers[4] = "all3d";
  
  // =====================================================================
  // configuration parameters

  std::string outfname   = pset.get<std::string>("OutputAnaFile");

 // =====================================================================

  
  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* evtree = new TTree("evflashmatch", "Event Flash-Match analysis");
  TTree* trtree = new TTree("trackflashmatch", "Track Flash-Match analysis");  
  //TTree* mcxingpt_tree = new TTree("mcxingptana", "Info on MC Crossing Point");
  //TTree* mcxingpt_prefilter_tree  = new TTree("mcxingptana_prefilter", "Info on MC Crossing Point");      

  // Event Index
  int run, subrun, event;
  bool ismc = true;
  bool badch_in_file = true;
  std::string inputimgs = "wire";
  int tag_neighborhood = 1;
  int verbosity = 0;
  std::string trigname = "triggersim";
  float fthreshold = 10;
  float fvertex_radius = 3;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;

  // Saved per event
  // ---------------
  // Truth Pixel Quantities
  truthdata.bindToTree( evtree );
  int nnu_pixels;        // number of neutrino pixels by most-nu-like cluster
  float fracnu_pixels;   // fraction of neutrino pixels covered by most-nu-like cluster
  int vertex_in_croi;       // found a good vertex
  int ngoodcroi;
  int nbadcroi;
  int num_rois;
  evtree->Branch( "nnu_pixels",     &nnu_pixels,     "nnu_pixels/I" );
  evtree->Branch( "fracnu_pixels",  &fracnu_pixels,  "fracnu_pixels/F" );
  evtree->Branch( "vertex_in_croi", &vertex_in_croi, "vertex_in_croi/I" );
  evtree->Branch( "ngoodcroi",      &ngoodcroi,      "ngoodcroi/I" );
  evtree->Branch( "nbadcroi",       &nbadcroi,       "nbadcroi/I" );
  evtree->Branch( "num_rois",       &num_rois,       "num_rois/I" );

  // Saved per track
  // ---------------
  // Flash Info
  truthdata.bindToTree( trtree );
  int isnutrack;            // flag indicating track is nu (1) or cosmic (0)
  float peratio;            // ratio of data/hypothesis
  float dcosmicchi2;        // Chi2(intime) - Chi2(cosmic)
  float chi2_intime;        // Chi2(intime)
  float chi2_cosmic;        // chi2(cosmic)
  float xmin;               // xmin of track (pre-extension from larlite)
  float xmax;               // xmax of track (pre-extension from larlite)
  trtree->Branch( "isnutrack",   &isnutrack,   "isnutrack/I" );
  trtree->Branch( "peratio",     &peratio,     "peratio/F" );
  trtree->Branch( "dcosmicchi2", &dcosmicchi2, "dcosmicchi2/F" );  
  trtree->Branch( "chi2_intime", &chi2_intime, "chi2_intime/F" );
  trtree->Branch( "xmin",        &xmin,        "xmin/F" );
  trtree->Branch( "xmax",        &xmax,        "xmax/F" );

  const float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5;
  
  // ==================================================================================================
  // ALGORITHMS
  
  // Space Charge Corrections
  larlitecv::SpaceChargeMicroBooNE sce;

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

    //if ( ientry%10==0 || verbosity>0 ) {
    std::cout << "===================================================" << std::endl;
    std::cout << "Entry " << ientry << ": (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    std::cout << "===================================================" << std::endl;

    // initialize the output variables
    truthdata.clear();
    nnu_pixels    = 0;
    fracnu_pixels = 0;      

    // =======================================================================
    // RETRIEVE DATA
    
    // images
    // --------------
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
      // // or we have to make the badch image from a ChStatus object
      // larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "tpc" );
      // std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      // ev_badch = new larcv::EventImage2D;
      // ev_badch->Emplace( std::move(chstatus_img_v) );
    }
    const std::vector<larcv::Image2D>& imgs_v = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& segs_v = ev_segs->Image2DArray();    

    // selected CROI
    // -------------
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
    
    // get truth
    // ---------
    larlite::event_mctruth* ev_mctruth   = NULL;
    larlite::event_mctrack* ev_mctrack   = NULL;
    larlite::event_mcshower* ev_mcshower = NULL;
    larlite::trigger* ev_trigger         = (larlite::trigger*) dataco[kSource].get_larlite_data(larlite::data::kTrigger, trigname);
    if ( ismc ) {
      ev_mctruth = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
      ev_mctrack = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
      ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");
    }

    // get tagger pixels/tracks
    // -------------------------
    larcv::EventPixel2D* ev_pix[kNumStages];
    larlite::event_track* ev_track[kNumStages];
    for (int i=0; i<kNumStages; i++) {
      try {
	ev_pix[i]    = (larcv::EventPixel2D*) dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D, stages_pixel_producers[i]);
	ev_track[i]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[i]);      
      }
      catch (...) {
	ev_pix[i] = NULL;
	ev_track[i] = NULL;
      }
    }

    // get user data, which we filled with flashmatch variables
    larlite::event_user* ev_flashdata = (larlite::event_user*)dataco[kCROIfile].get_larlite_data(larlite::data::kUserInfo,"croicutresults");
    larlite::user_info& flashdata = ev_flashdata->at(0);
    
    // =======================================================================
    // MC INFO ANALYSIS

    // extract the truth quantities of interest
    if ( ismc ) {
      larlitecv::extractTruth( *ev_mctruth, *ev_mctrack, truthdata );
    }

    // quantities filled if MC present
    std::vector<larcv::Image2D> nupix_imgs_v;
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
      std::cout << "Vertex 3D Coordinates (uncorrected):      (" << dpos[0] << "," << dpos[1] << "," << dpos[2] << ")" << std::endl;

      // did any of the ROIs contain the vertex?
      vertex_in_croi = 0;
      ngoodcroi = 0;
      nbadcroi  = 0;
      for ( auto& roi : containedrois_v ) {
	int nplanes_in_roi = 0;
	for (size_t p=0; p<3; p++ ) {
	  const larcv::ImageMeta& bb = roi.BB( (larcv::PlaneID_t)p );
	  if ( vertex_tick>=bb.min_y() && vertex_tick<=bb.max_y() && vertex_col[p]>=bb.min_x() && vertex_col[p]<=bb.max_x() )
	    nplanes_in_roi++;
	}
	if (nplanes_in_roi>=2) {
	  vertex_in_croi = 1;
	  ngoodcroi++;
	}
	else {
	  nbadcroi++;
	}
      }
    }//if is mc
      
    // make truth pixel counts
    // -----------------------
    // count the pixels. determine if cosmic and neutrino are tagged. also if neutrino is in rois
    // we loop through the rows and cols
    int num_nu_pixels_tot = 0;
    if ( ismc ) {
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
	      //nvertex_pixels[p]++;
	      //nvertex_pixels[3]++;
	    }
	      
	    // above threshold. is it a neutrino pixel?
	    bool in_seg_image = false;
	    int seg_row = -1;
	    int seg_col = -1;
	      
	    const larcv::Image2D& segimg = segs_v.at(p);
	    float x = imgs_v.at(p).meta().pos_x(col);
	    float y = imgs_v.at(p).meta().pos_y(row);
	    if ( x>segs_v.at(p).meta().min_x() && x<segs_v.at(p).meta().max_x()
		 && y>segs_v.at(p).meta().min_y() && y<segs_v.at(p).meta().max_y() ) {
	      in_seg_image = true;
	      seg_row = segs_v.at(p).meta().row(y);
	      seg_col = segs_v.at(p).meta().col(x);
	    }

	    if ( in_seg_image && segs_v.at(p).pixel(seg_row,seg_col)>0 ) {		
	      is_nu_pix = true;
	      num_nu_pixels_tot++;
	    }
	    
	    if ( is_nu_pix )
	      nupix_img.set_pixel( row, col, 1.0 );
	    if ( near_vertex )
	      nupix_img.set_pixel( row, col, 10.0 );
	    
	  }//end of col loop
	}//end of row loop
	nupix_imgs_v.emplace_back( std::move(nupix_img) );
      }//end of loop over planes for counting neutrino/cosmic pixels
    }//end of MC Pixel Counts
    

    // loop over tagger stage results
    // we count the number neutrino pixels a tagger track overlaps
    // we use the largest one to tag the "neutrino" cluster
    int istage=kAll;
    if ( ev_pix[istage]==NULL )
      throw std::runtime_error("All stage pixels/tracks not saved");
    
    int numtracks = ev_pix[istage]->Pixel2DClusterArray(0).size();
    std::vector<int> num_tagged_nupixels(numtracks,0);
      
    for ( size_t itrack=0; itrack<numtracks; itrack++ ) {
      
      for ( size_t p=0; p<3; p++ ) {
	
	auto const& pixcluster  = ev_pix[istage]->Pixel2DClusterArray(p).at(itrack);
	
	// for each track cluster, we loop over pixels and store unique pixels by using a set
	std::set< std::vector<int> > tagged_set;
	for ( auto const& pix : pixcluster ) {
	  
	  bool tagged = false;
	  for ( int dr=-tag_neighborhood; dr<=tag_neighborhood; dr++ ) {
	    if ( tagged ) break;
	    int r = pix.Y()+dr;
	    if ( r<0 || r>=(int)nupix_imgs_v.at(p).meta().rows() ) continue;
	    for ( int dc=-tag_neighborhood; dc<=tag_neighborhood; dc++ ) {
	      if ( tagged ) break;
	      int c = pix.X()+dc;
	      if ( c<0 || c>=(int)nupix_imgs_v.at(p).meta().cols() ) continue;
	      if ( nupix_imgs_v.at(p).pixel(r,c)>0.0 )
		tagged = true;
	    }//dc loop
	  }//dr loop
	  if ( tagged )
	    num_tagged_nupixels[itrack]++;
	}//end of pixel loop
      }// number of planes
    }// number of tracks

    // get largest neutrino track
    int max_itrack  = -1;
    int max_ntagged = 0;
    for ( size_t itrack=0; itrack<numtracks; itrack++ ) {
      if ( max_itrack<num_tagged_nupixels[itrack] ) {
	max_itrack  = itrack;
	max_ntagged = num_tagged_nupixels[itrack];
      }
    }

    if ( max_itrack<0 ) {
      nnu_pixels = -1;
      fracnu_pixels = -1;
    }
    else {
      nnu_pixels = max_ntagged;
      fracnu_pixels = float(nnu_pixels)/float(num_nu_pixels_tot);
    }

    // correlate flash-match variables with (neutrino/cosmic) tag
    std::vector<double>* flashdata_chi2        = flashdata.get_darray("flashmatch_intimechi2");
    std::vector<double>* flashdata_cosmicratio = flashdata.get_darray("cosmicratio_dchi2");
    std::vector<double>* flashdata_peratio     = flashdata.get_darray("totalpe_peratio");      
    std::vector<double>* flashdata_dwall       = flashdata.get_darray("containment_dwall");

    if ( !flashdata_chi2 || !flashdata_cosmicratio || !flashdata_peratio || !flashdata_dwall ) {
      std::cout << "Empty result" << std::endl;
      evtree->Fill();
      continue;
    }
      
    for ( size_t itrack=0; itrack<numtracks; itrack++ ) {
      // tag if neutrino or cosmic
      if ( itrack==max_itrack )
	isnutrack = 1;
      else
	isnutrack = 0;
      // copy over variables
      peratio = flashdata_peratio->at(itrack);
      dcosmicchi2 = flashdata_cosmicratio->at(itrack);
      if ( std::isinf(dcosmicchi2) )
	dcosmicchi2 = 1e6;
      
      std::cout << "dcosmicchi2=" << dcosmicchi2 << std::endl;
      chi2_intime = flashdata_chi2->at(itrack);
      const larlite::track& lltrack = ev_track[kAll]->at(itrack);
      xmin = 1e9;
      xmax = 0;
      for ( int ipt=0; ipt<lltrack.NumberTrajectoryPoints(); ipt++ ) {
	float x = lltrack.LocationAtPoint(ipt).X();
	if ( xmin>x )
	  xmin = x;
	if ( xmax<x )
	  xmax = x;
      }
      trtree->Fill();
    }

    //================================================================================

    
    evtree->Fill();
    
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
