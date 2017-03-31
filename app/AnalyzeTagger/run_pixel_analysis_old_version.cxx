#include "Base/DataCoordinator.h"

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
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "dwall.h"
#include "extractTruthMethods.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "GapChs/EmptyChannelAlgo.h"


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

  enum { kSource=0, kCROIfile, kNumSourceTypes } SourceTypes_t;
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "AnalyzeTagger" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }

  // configuration parameters
  float fthreshold       = pset.get<float>("PixelThreshold");
  int tag_neighborhood   = pset.get<int>("TagNeighborhood",10);
  int fvertex_radius     = pset.get<int>("PixelRadius");
  int verbosity          = pset.get<int>("Verbosity",0);

  // setup output
  enum { kThruMu=0, kStopMu, kUntagged, kCROI, kNumStages } Stages_t;


  TFile* rfile = new TFile("output_pixel_analysis_rerunpad20.root", "recreate");
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

  // Crossing Point tags
  int true_intime_thrumu;
  int true_intime_stopmu;
  int true_crossingpoints;
  int tagged_crossingpoints;
  int proposed_crossingpoints;

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
  tree->Branch("vtx_in_croi",             &vertex_in_croi,         "vtx_in_croi/I" );
  tree->Branch("true_intime_thrumu",      &true_intime_thrumu,     "true_intime_thrumu/I");
  tree->Branch("true_intime_stopmu",      &true_intime_stopmu,     "true_intime_stopmu/I");
  tree->Branch("true_crossingpoints",     &true_crossingpoints,     "true_crossingpoints/I");
  tree->Branch("tagged_crossingpoints",   &tagged_crossingpoints,   "tagged_crossingpoints/I");
  tree->Branch("proposed_crossingpoints", &proposed_crossingpoints, "proposed_crossingpoints/I");

  // Space Charge Corrections
  larlitecv::SpaceChargeMicroBooNE sce;

  // Empty Channel Algo
  larlitecv::EmptyChannelAlgo emptyalgo;

  int nentries = dataco[kCROIfile].get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(run,subrun,event);

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
    proposed_crossingpoints = 0;
    tagged_crossingpoints = 0;
    true_crossingpoints = 0;
    true_intime_stopmu = 0;
    true_intime_thrumu = 0;
    vertex_in_croi = 0;


    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_segs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"modimg");
    larcv::EventImage2D* ev_badch  = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"gapchs");

    larcv::EventPixel2D* ev_pix[kNumStages] = {0};
    ev_pix[kThruMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"thrumupixels");
    ev_pix[kStopMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"stopmupixels");
    ev_pix[kUntagged]  = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"untaggedpixels");
    ev_pix[kCROI]      = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"croipixels");

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>& segs_v   = ev_segs->Image2DArray();

    // make badch image from ChStatus object
    larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "tpc" );
    std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );

    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    const std::vector<larcv::ROI>& containedrois_v = ev_contained_roi->ROIArray();
    num_rois = (int)containedrois_v.size();
    std::cout << "====ROIs===========================" << std::endl;
    for ( auto& roi : containedrois_v ) {
      std::cout << " roi: " << roi.dump();
    }
    std::cout << "===================================" << std::endl;

    // get the boundary end point info
    larcv::EventPixel2D* ev_spacepoints[7];
    std::string spacepoint_producers[7] = { "topspacepts", "botspacepts", "upspacepts", "downspacepts", "anodepts", "cathodepts", "imgendpts" };
    for ( int i=0; i<7; i++ ) {
      ev_spacepoints[i] = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,spacepoint_producers[i]);
    }

    // get other information, e.g. truth
    larlite::event_mctruth* ev_mctruth   = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
    larlite::event_mctrack* ev_mctrack   = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");
    larlite::trigger* ev_trigger         = (larlite::trigger*)       dataco[kSource].get_larlite_data(larlite::data::kTrigger,"triggersim");

    // extract the truth quantities of interest
    larlitecv::extractTruth( *ev_mctruth, *ev_mctrack, truthdata );

    // get the vertex in the pixel coordinates
    std::vector<int> wid(3,-1);
    std::vector<int> vertex_col(3,-1);
    std::vector<double> dpos(3);
    for (int i=0; i<3; i++ ) dpos[i] = truthdata.pos[i];
    std::vector<double> vtx_offset = sce.GetPosOffsets( dpos[0], dpos[1], dpos[2] );
    std::vector<double> vtx_sce(3,0);
    for (int i=1; i<3; i++ ) vtx_sce[i] = dpos[i] + vtx_offset[i];
    vtx_sce[0] = dpos[0] - vtx_offset[0] + 0.7;
    
    for (size_t p=0; p<3; p++) {
      wid[p] = ::larutil::Geometry::GetME()->WireCoordinate( vtx_sce, p );
      if ( wid[p]>=0 && wid[p]<3456 )
	vertex_col[p] = imgs_v.at(p).meta().col(wid[p]);
    }
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float vertex_tick = vtx_sce[0]/cm_per_tick + 3200.0;
    int vertex_row = -1;
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
      if (nplanes_in_roi==3) {
	vertex_in_croi = 1;
	break;
      }
    }

    // loop over MC tracks, get end points of muons
    std::vector< std::vector<int> > start_pixels;
    std::vector< std::vector<float> > start_crossingpts;
    std::vector< std::vector<int> > end_pixels;
    std::vector< std::vector<float> > end_crossingpts;
    int intime_cosmics = 0;
    for ( auto const& track : *ev_mctrack ) {
      if ( std::abs(track.PdgCode())!=13  ) continue;
      if ( track.size()==0 ) continue;
      const TLorentzVector& track_start = track.front().Position();
      std::vector<float> fstart(3,0);
      fstart[0] = track_start.X();
      fstart[1] = track_start.Y();
      fstart[2] = track_start.Z();
      std::vector<float> fend(3,0);
      fend[0]   = track.back().Position().X();
      fend[1]   = track.back().Position().Y();
      fend[2]   = track.back().Position().Z();

      int track_start_boundary = 0;
      float track_start_dwall = larlitecv::dwall( fstart, track_start_boundary );
      int track_end_boundary = 0;
      float track_end_dwall = larlitecv::dwall( fend, track_end_boundary );

      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();

      // space charge corrections
      std::vector<double> start_offset = sce.GetPosOffsets( first_step.X(), first_step.Y(), first_step.Z() );
      std::vector<double> end_offset   = sce.GetPosOffsets( last_step.X(),  last_step.Y(),  last_step.Z() );
      //std::vector<double> start_offset(3,0);
      //std::vector<double> end_offset(3,0);

      std::vector<float> sce_start(3);
      sce_start[0] = first_step.X()-start_offset[0]+0.7;
      sce_start[1] = first_step.Y()+start_offset[1];
      sce_start[2] = first_step.Z()+start_offset[2];
      Double_t sce_start_xyz[3] = { sce_start[0], sce_start[1], sce_start[2] };

      std::vector<float> sce_end(3);
      sce_end[0] = last_step.X()-end_offset[0]+0.7;
      sce_end[1] = last_step.Y()+end_offset[1];
      sce_end[2] = last_step.Z()+end_offset[2];
      Double_t sce_end_xyz[3] = { sce_end[0], sce_end[1], sce_end[2] };

      std::vector< int > start_pix(4); // (row, u-plane, v-plane, y-plane)
      std::vector< int > end_pix(4); // (row, u-plane, v-plane, y-plane)
      start_pix[0] = (first_step.T()*1.0e-3 - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_start[0]/cm_per_tick + 3200;
      end_pix[0]   = (last_step.T()*1.0e-3  - (ev_trigger->TriggerTime()-4050.0))/0.5 + sce_end[0]/cm_per_tick   + 3200;
      for ( size_t p=0; p<3; p++ ) {
        start_pix[p+1] = larutil::Geometry::GetME()->WireCoordinate( sce_start_xyz, p );
        end_pix[p+1]   = larutil::Geometry::GetME()->WireCoordinate( sce_end_xyz,   p );
      }

      bool start_intime = false;
      bool start_crosses = false;
      if ( start_pix[0]>imgs_v.front().meta().min_y() && start_pix[0]<imgs_v.front().meta().max_y() ) {
        start_intime = true;
        start_pix[0] = imgs_v.front().meta().row( start_pix[0] );

        if ( track_start_dwall < 10.0 ) {
          start_pixels.emplace_back( std::move(start_pix) );
          start_crossingpts.emplace_back( std::move(sce_start) );
          start_crosses = true;
          true_crossingpoints++;
        }
      }

      bool end_intime = false;
      bool end_crosses = false;
      if ( end_pix[0]>imgs_v.front().meta().min_y() && end_pix[0]<imgs_v.front().meta().max_y() ) {
        end_intime = true;
        end_pix[0]   = imgs_v.front().meta().row( end_pix[0] );
        if ( track_end_dwall < 10.0 ) {
          end_pixels.emplace_back( std::move(end_pix) );
          end_crossingpts.emplace_back( std::move(sce_end) );
          end_crosses = true;
          true_crossingpoints++;
        }
      }

      if ( start_intime || end_intime ) {
        intime_cosmics++;
        if ( start_intime && !start_crosses ) {
          std::cout << "start point: (" << fstart[0] << "," << fstart[1] << "," << fstart[2] << ")"
                    << " dwall=" << track_start_dwall << " intime=" << start_intime
                    << " tick=" << imgs_v.front().meta().pos_y( start_pix[0] )
                    << " row=" << start_pix[0]
                    << std::endl;
          //throw std::runtime_error("start point does not cross boundary?");
        }
        else if ( start_crosses && end_crosses )
          true_intime_thrumu++;
        else if ( start_crosses && !end_crosses )
          true_intime_stopmu++;
      }

    }//end of loop over mc tracks

    std::cout << "number of intime cosmics: "       << intime_cosmics << std::endl;
    std::cout << "number of intime thrumu: "        << true_intime_thrumu << std::endl;
    std::cout << "number of intime stopmu: "        << true_intime_stopmu << std::endl;
    std::cout << "number of true crossing points: " << true_crossingpoints << std::endl;

    // make truth pixel counts
    // count the pixels. determine if cosmic and neutrino are tagged. also if neutrino is in rois
    // we loop through the rows and cols
    std::vector<larcv::Image2D> nupix_imgs_v;
    for (size_t p=0; p<3; p++) {

      // we create a neutrino pixel image, to make things easier downstream
      larcv::Image2D nupix_img( imgs_v.at(p).meta() );
      nupix_img.paint(0);
      for (size_t row=0; row<imgs_v.at(p).meta().rows(); row++) {
        for (size_t col=0; col<imgs_v.at(p).meta().cols(); col++) {

	        // check if this is a pixel of interest
          if ( imgs_v.at(p).pixel(row,col)<fthreshold ) continue;

	        nthreshold_pixels[p]++;
	        nthreshold_pixels[3]++;

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
          const larcv::Image2D& segimg = segs_v.at(p);
          float x = imgs_v.at(p).meta().pos_x(col);
          float y = imgs_v.at(p).meta().pos_y(row);
          bool in_seg_image = false;
          int seg_row = -1;
          int seg_col = -1;
          if ( x>segs_v.at(p).meta().min_x() && x<segs_v.at(p).meta().max_x()
            && y>segs_v.at(p).meta().min_y() && y<segs_v.at(p).meta().max_y() ) {
            in_seg_image = true;
            seg_row = segs_v.at(p).meta().row(y);
            seg_col = segs_v.at(p).meta().col(x);
          }
          if ( in_seg_image && segs_v.at(p).pixel(seg_row,seg_col)>0 ) {

            is_nu_pix = true;
            nnu_pixels[p]++;
            nnu_pixels[3]++;
            //if (near_vertex) {
            //  nvertex_pixels[p]++;
            //  nvertex_pixels[3]++;
	          //}

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
              if ( r<0 || r>=imgs_v.at(p).meta().rows() ) continue;
              for ( int dc=-tag_neighborhood; dc<=tag_neighborhood; dc++ ) {
                int c = pix.X()+dc;
                if ( c<0 || c>=imgs_v.at(p).meta().cols() ) continue;
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
            bool isnupix  = ( nupix_imgs_v.at(p).pixel(r,c)>=1.0  ) ? true : false;
            bool isvertex = ( nupix_imgs_v.at(p).pixel(r,c)>=10.0 ) ? true : false;

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

    // did any of the ROIs catch neutrino pixels?
    for ( int r=vertex_row-fvertex_radius; r<=vertex_row+fvertex_radius; r++ ) {
      if ( r<0 || r>=imgs_v.front().meta().rows() ) continue;
      for (size_t p=0; p<3; p++) {
        const larcv::ImageMeta& meta = imgs_v.at(p).meta();
        for ( int c=vertex_col[p]-fvertex_radius; c<=vertex_col[p]+fvertex_radius; c++ ) {
          if ( c<0 || c>=meta.cols() ) continue;

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

    // analyze proposed boundary points
    std::cout << "Analyze Boundary Points" << std::endl;
    proposed_crossingpoints = 0;
    for (int i=0; i<7; i++) {
      if ( ev_spacepoints[i]==NULL)
        throw std::runtime_error("wtf");
      std::cout << " endtype " << spacepoint_producers[i] << ": ";
      std::cout << ev_spacepoints[i]->Pixel2DArray(0).size() << std::endl;
      proposed_crossingpoints += (int)(ev_spacepoints[i]->Pixel2DArray(0).size());
    }


    std::cout << "Match Truth to Tagged" << std::endl;

    std::vector<bool> matched_startpoint( start_pixels.size(), false );
    std::vector<bool> matched_endpoint( end_pixels.size(), false );

    std::vector< std::vector<int> >* p_crossing_pixel_v[2] = { &start_pixels, &end_pixels }; // truth pixels for crossing points
    std::vector< std::vector<float> >* p_crossingpts_v[2]  = { &start_crossingpts, &end_crossingpts }; // truth 3D positions
    std::vector<bool>* p_matched_v[2] = { &matched_startpoint, &matched_endpoint };


    for ( int v=0; v<2; v++ ) {
      for ( int ipix=0; ipix<(int)p_crossing_pixel_v[v]->size(); ipix++ ) {

        // we need to get the 3D position to compare against
        const std::vector<int>&  pixinfo     = p_crossing_pixel_v[v]->at(ipix); // position in image
        const std::vector<float>& crossingpt = p_crossingpts_v[v]->at(ipix);    // position in 3D

        // use TPC position to get X
        std::vector<float> crossingpt_tpcx(3);
        crossingpt_tpcx[0] = (imgs_v.at(0).meta().pos_y( pixinfo[0] )-3200.0)*cm_per_tick;
        crossingpt_tpcx[1] = crossingpt[1];
        crossingpt_tpcx[2] = crossingpt[2];

        // scan for pixel, loop over types and pts
        bool matched = false;
        for (int i=0; i<7; i++) {
          if ( matched )
            break;

          for ( int j=0; j<(int)ev_spacepoints[i]->Pixel2DArray(0).size(); j++ ) {


            std::vector<float> intersect(2,0.0);
            std::vector<int> wids(3,0);
            int crossing = 0;
            double triangle_area = 0.0;
            for (int p=0; p<3; p++) {
              wids[p] = ev_spacepoints[i]->Pixel2DArray(p).at(j).X();
            }


            larcv::UBWireTool::wireIntersection( wids, intersect, triangle_area, crossing );
            float x = ( imgs_v.at(0).meta().pos_y( ev_spacepoints[i]->Pixel2DArray(0).at(j).Y() ) - 3200.0 )*cm_per_tick;

            std::vector<float> spacepoints(3);
            spacepoints[0] = x;
            spacepoints[1] = intersect[1];
            spacepoints[2] = intersect[0];

            float dist = 0;
            for (int d=0; d<3; d++) {
              dist += (spacepoints[d]-crossingpt_tpcx[d])*(spacepoints[d]-crossingpt_tpcx[d]);
            }
            dist = sqrt(dist);
            //std::cout << "true[" << v << "," << ipix << "] vs. proposed[" << i << "," << j << "] dist=" << dist << std::endl;
            if (dist<20.0) {
              matched = true;
            }

            if ( matched )
              break;
          }// end of loop over tagged points of type i
        }//end of boundary point types

        p_matched_v[v]->at(ipix) = matched;

        if ( matched )
          tagged_crossingpoints++;
      }
    }

    std::cout << "Proposed Crossing Points: " << proposed_crossingpoints << std::endl;
    std::cout << "Tagged Crossing Points: " << tagged_crossingpoints << std::endl;
    std::cout << "True Crossing Points: " << true_crossingpoints << std::endl;

#ifdef USE_OPENCV
    // draw image
    for ( size_t p=0; p<cvleftover_v.size(); p++ ) {
      auto& leftover = cvleftover_v.at(p);

      // draw truth end points!
      for ( auto const& start_pix : start_pixels) {
        cv::circle( leftover, cv::Point(start_pix[p+1],start_pix[0]), 4, cv::Scalar( 0, 255, 0 ), 2, -1 );
      }
      for ( auto const& end_pix : end_pixels) {
        cv::circle( leftover, cv::Point(end_pix[p+1],end_pix[0]), 4, cv::Scalar( 255, 255, 0 ), 2, -1 );
      }

      // draw proposed end points
      for ( int i=0; i<7; i++) {
        for ( auto const& endpt : ev_spacepoints[i]->Pixel2DArray(p) ) {
          cv::drawMarker( leftover, cv::Point(endpt.X(), endpt.Y()),  cv::Scalar(0, 0, 255), cv::MARKER_CROSS, 6, 2);
        }
      }

      // draw roi
      for ( auto const& roi : containedrois_v ) {
        larcv::draw_bb( leftover, imgs_v.front().meta(), roi.BB(p), 255, 0, 255, 2 );
      }

      // sce vertex
      cv::circle( leftover, cv::Point(vertex_col[p],vertex_row), 4, cv::Scalar(0,0,255),   2, -1 );
      cv::circle( leftover, cv::Point(vertex_col[p],vertex_row), 3, cv::Scalar(0,255,255), 1, -1 );      

      std::stringstream ss;
      ss << "leftover_i" << ientry << "_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
      std::cout << "write: " << ss.str() << std::endl;
      cv::imwrite( ss.str(), leftover );
    }
#endif

    tree->Fill();

  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
