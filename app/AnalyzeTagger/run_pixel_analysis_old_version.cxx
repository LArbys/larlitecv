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
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

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


// OpenCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"

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


  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
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

  // Tagged Pixel Quantities
  int ncosmic_tagged[kNumStages][4];       // number of non-neutrino pixels tagged
  int ncosmic_tagged_once[kNumStages][4];  // number of non-neutrino pixels tagged
  int ncosmic_tagged_many[kNumStages][4];  // number of non-neutrino pixels tagged
  int nnu_tagged[kNumStages][4];           // number of neutrino pixels tagged
  int nvertex_tagged[kNumStages][4];       // number of neutrino pixels within some pixel radius of vertex tagged
  std::stringstream s_arr;
  s_arr << "[" << (int)kNumStages << "][4]/I";

  // ROI quantities
  int num_rois;     // number of identified ROis
  int nnu_inroi[4]; // number of nu pixels contained in the CROI

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

  tree->Branch("ncosmic_tagged",      ncosmic_tagged,      std::string("ncosmic_tagged"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_once", ncosmic_tagged_once, std::string("ncosmic_tagged_once"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_many", ncosmic_tagged_many, std::string("ncosmic_tagged_many"+s_arr.str()).c_str() );
  tree->Branch("nnu_tagged",          nnu_tagged,          std::string("nnu_tagged"+s_arr.str()).c_str() );
  tree->Branch("nvertex_tagged",      nvertex_tagged,      std::string("nvertex_tagged"+s_arr.str()).c_str() );

  tree->Branch("num_rois", &num_rois, "num_rois/I");
  tree->Branch("nnu_inroi", nnu_inroi, "nnu_inroi[4]" );

  // Space Charge Corrections
  larlitecv::SpaceChargeMicroBooNE sce;

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
      for (int istage=0; istage<kNumStages; istage++) {
	      ncosmic_tagged[istage][p] = 0;
	      ncosmic_tagged_once[istage][p] = 0;
	      ncosmic_tagged_many[istage][p] = 0;
	      nnu_tagged[istage][p] = 0;
	      nvertex_tagged[istage][p] = 0;
      }
    }

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

    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    const std::vector<larcv::ROI>& containedrois_v = ev_contained_roi->ROIArray();
    num_rois = (int)containedrois_v.size();
    std::cout << "====ROIs===========================" << std::endl;
    for ( auto& roi : containedrois_v ) {
      std::cout << " roi: " << roi.dump();
    }
    std::cout << "===================================" << std::endl;

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
    for (size_t p=0; p<3; p++) {
      wid[p] = ::larutil::Geometry::GetME()->WireCoordinate( dpos, p );
      if ( wid[p]>=0 && wid[p]<3456 )
	      vertex_col[p] = imgs_v.at(p).meta().col(wid[p]);
    }
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float ftick = truthdata.pos[0]/cm_per_tick + 3200.0;
    int vertex_row = -1;
    if ( ftick >= imgs_v.at(0).meta().min_y() && ftick<=imgs_v.at(0).meta().max_y() )
      vertex_row = imgs_v.at(0).meta().row( ftick );


    std::cout << "Vertex Pixel Coordinates (" << vertex_row << ", " << vertex_col[0] << "," << vertex_col[1] << "," << vertex_col[2] << ")" << std::endl;

    // loop over MC tracks, get end points of muons
    std::vector< std::vector<int> > start_pixels;
    std::vector< std::vector<int> > end_pixels;
    int intime_cosmics = 0;
    for ( auto const& track : *ev_mctrack ) {
      if ( std::abs(track.PdgCode())!=13  ) continue;
      if ( track.size()==0 ) continue;
      const TLorentzVector& track_start = track.front().Position();
      std::vector<float> fstart(3);
      fstart[0] = track_start.X();
      fstart[1] = track_start.Y();
      fstart[2] = track_start.Z();

      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();

      // space charge corrections
      std::vector<double> start_offset = sce.GetPosOffsets( first_step.X(), first_step.Y(), first_step.Z() );
      std::vector<double> end_offset   = sce.GetPosOffsets( last_step.X(),  last_step.Y(),  last_step.Z() );

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

      if ( start_pix[0]>imgs_v.front().meta().min_y() && start_pix[0]<imgs_v.front().meta().max_y() ) {
        start_pix[0] = imgs_v.front().meta().row( start_pix[0] );
        start_pixels.emplace_back( std::move(start_pix) );
        intime_cosmics++;
      }

      if ( end_pix[0]>imgs_v.front().meta().min_y() && end_pix[0]<imgs_v.front().meta().max_y() ) {
        end_pix[0]   = imgs_v.front().meta().row( end_pix[0] );
        end_pixels.emplace_back( std::move(end_pix) );
      }

    }//end of loop over mc tracks

    std::cout << "number of intime cosmics: "<< intime_cosmics << std::endl;


    // make thruth pixel counts
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
    std::vector<larcv::Image2D> leftover_v;
    std::vector<cv::Mat> cvleftover_v;
    for ( auto const& img : imgs_v ) {
      larcv::Image2D leftover(img);
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( leftover, fthreshold, 100.0 );
      leftover_v.emplace_back(leftover);
      cvleftover_v.emplace_back(cvimg);
    }

    for ( int istage=0; istage<kNumStages; istage++ ) {
      for ( size_t p=0; p<3; p++ ) {
        // we need images to track the number of times pixels are Tagged
        larcv::Image2D img_counter( imgs_v.at(p).meta() );
        img_counter.paint(0);

        for ( auto const& pixcluster : ev_pix[istage]->Pixel2DClusterArray( p ) ) {

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
          for ( auto const& pix : tagged_set ) {
            float val = img_counter.pixel( pix[1], pix[0]) + 1.0;
            img_counter.set_pixel( pix[1], pix[0], val );
            leftover_v.at(p).set_pixel(pix[1],pix[0],0);
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
            else {
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[0] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[1] = 0;
              cvleftover_v.at(p).at<cv::Vec3b>(cv::Point(pix[0],pix[1]))[2] = 200;
            }
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
          }
        }
      }
    }

    // draw image
    int p = 0;
    for ( auto const& leftover : cvleftover_v ) {

      // draw truth end points!
      for ( auto const& start_pix : start_pixels) {
        cv::circle( leftover, cv::Point(start_pix[p+1],start_pix[0]), 4, cv::Scalar( 0, 255, 0 ), 2, -1 );
      }
      for ( auto const& end_pix : end_pixels) {
        cv::circle( leftover, cv::Point(end_pix[p+1],end_pix[0]), 4, cv::Scalar( 255, 255, 0 ), 2, -1 );
      }

      std::stringstream ss;
      ss << "leftover_i" << ientry << "_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
      cv::imwrite( ss.str(), leftover );
      p++;
    }


    tree->Fill();

    break;

  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
