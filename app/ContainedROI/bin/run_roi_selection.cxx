#include <iostream>
#include <string>
#include <vector>

// larcv
#include "DataFormat/EventROI.h"
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventImage2D.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/opflash.h"
#include "DataFormat/chstatus.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "ThruMu/EmptyChannelAlgo.h"
#include "ContainedROIConfig.h"
#include "ContainedROI.h"
#include "FlashROIMatching.h"

#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif


int main( int nargs, char** argv ) {

  std::cout << "[[ ROI SELECTION ]]" << std::endl;

  // we need to have a data coordinator for each stage because the number of entries could be different.
  // we'll coordinate by using event,subrun,run information
  std::string data_folder = "~/working/data/larbys/cosmic_tagger_dev/";

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
  dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv.root",   "larcv" );   //stopmu-tagger output
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larlite.root", "larlite" ); //stopmu-tagger output

  larlitecv::DataCoordinator dataco_output;
  dataco_output.set_outputfile( "output_containedroi_larcv.root", "larcv");
  dataco_output.set_outputfile( "output_containedroi_larlite.root", "larlite");

  dataco_source.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_thrumu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_stopmu.configure( "containedroi.cfg", "StorageManager", "IOManager", "ContainedROI" );
  dataco_output.configure( "containedroi.cfg", "StorageManagerOutput", "IOManagerOutput", "ContainedROI" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();
  dataco_output.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "containedroi.cfg" );
  larcv::PSet pset = cfg.get<larcv::PSet>("ContainedROI");
  larcv::PSet contained_pset = pset.get<larcv::PSet>("ContainedROIConfig");
  larcv::PSet flash_pset     = pset.get<larcv::PSet>("FlashROIMatchingConfig");
  std::vector<std::string> flashproducers = pset.get< std::vector<std::string> >("OpFlashProducers");

  larlitecv::ContainedROIConfig contained_cfg = larlitecv::CreateContainedROIConfig( contained_pset );
  larlitecv::FlashROIMatchingConfig flash_cfg = larlitecv::MakeFlashROIMatchingConfigFromFile( "containedroi.cfg"); // add pset interface

  // contained ROI selection algorithm
  larlitecv::ContainedROI algo( contained_cfg );
  larlitecv::FlashROIMatching flash_matching( flash_cfg );

  // event loop

  int nentries = dataco_stopmu.get_nentries("larcv");
  int user_nentries =   pset.get<int>("NumEntries",-1);
  int user_startentry = pset.get<int>("StartEntry",-1);
  int start_entry = 0;
  if ( user_startentry>=0 ) {
    start_entry = user_startentry;
  }
  int end_entry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+start_entry<nentries )
      end_entry = user_nentries+start_entry;
  }

  for ( int ientry=start_entry; ientry<end_entry; ientry++ ) {

    // load same event in all three data sets
    dataco_stopmu.goto_entry(ientry,"larcv");
    int run,subrun,event;
    dataco_stopmu.get_id(run,subrun,event);
    std::cout << "entry " << ientry << std::endl;
    std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    dataco_thrumu.goto_event(run,subrun,event,"larcv");
    dataco_source.goto_event(run,subrun,event,"larcv");
    dataco_output.set_id(run, subrun, event);

    algo.m_run = run;
    algo.m_subrun = subrun;
    algo.m_event = event;
    algo.m_entry = ientry;

    // get images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"modimgs");
    larcv::EventImage2D* ev_seg    = (larcv::EventImage2D*)dataco_source.get_larcv_data(larcv::kProductImage2D,"segment");
    larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"marked3d");
    larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D,"stopmu"); 
    larcv::EventROI* rois          = (larcv::EventROI*)    dataco_source.get_larcv_data(larcv::kProductROI, "tpc" );


    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
    const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();

    // use larlite chstatus to get badch
    larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco_source.get_larlite_data( larlite::data::kChStatus, "chstatus" );
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::cout << "ch status planes: " << ev_status->size() << std::endl;
    std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;

    // larlite opflash data
    std::vector< larlite::event_opflash* > opflash_containers;
    for ( auto &flashproducer : flashproducers ) {
      std::cout << "search for flash hits from " << flashproducer << ". ";
      larlite::event_opflash* opdata = (larlite::event_opflash*)dataco_source.get_larlite_data(larlite::data::kOpFlash, flashproducer );
      std::cout << " has " << opdata->size() << " flashes" << std::endl;
      opflash_containers.push_back( opdata );
    }

    // find untagged pixel clusters, matched in all three planes
    std::vector< std::vector<larcv::Pixel2DCluster> > untagged_pixel_clusters;
    algo.SetTruthROI( rois->ROIArray().at(0) );
    std::vector<larcv::ROI> untagged_rois = algo.SelectROIs( imgs_v, thrumu_v, stopmu_v, badch_v, untagged_pixel_clusters );

    // retrieve the thrumu and stopmu clusters as well
    larcv::EventPixel2D*  ev_thrumu_pixels = (larcv::EventPixel2D*) dataco_thrumu.get_larcv_data(larcv::kProductPixel2D,   "thrumu2d" );
    larlite::event_track* ev_thrumu_tracks = (larlite::event_track*)dataco_thrumu.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    larcv::EventPixel2D*  ev_stopmu_pixels = (larcv::EventPixel2D*) dataco_stopmu.get_larcv_data(larcv::kProductPixel2D,   "stopmupixels");
    larlite::event_track* ev_stopmu_tracks = (larlite::event_track*)dataco_stopmu.get_larlite_data( larlite::data::kTrack, "stopmutracks");

    // we pass these clusters to the t
    larcv::Image2D yplane_seg = ev_seg->Image2DArray().at(2); // just copy;
    flash_matching.setNeutrinoYPlaneSegmentImage( &yplane_seg );
    flash_matching.setRSE( ientry, run, subrun, event );
    std::vector<larcv::ROI> flash_matched_rois = flash_matching.SelectFlashConsistentROIs( opflash_containers, imgs_v, 
      untagged_pixel_clusters, untagged_rois, *ev_thrumu_pixels, *ev_stopmu_pixels );


    // visualize the output
    std::vector< cv::Mat > cvimgs_v;
    for (size_t p=0; p<3; p++) {
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( imgs_v.at(p), 5.0, 50.0 );

      // draw interaction ROI
      if (rois->ROIArray().size()>0 && p<(int)rois->ROIArray().at(0).BB().size() && p<(int)imgs_v.size() ) {
       larcv::draw_bb( cvimg, imgs_v.at(p).meta(), rois->ROIArray().at(0).BB().at(p), 200, 0, 200, 2 );
      }

      // label thrumu and stopmu pixels
      int nthrumu = ev_thrumu_pixels->Pixel2DClusterArray(p).size();
      for ( int imu=0; imu<nthrumu; imu++ ) {
        const larcv::Pixel2DCluster& cluster = ev_thrumu_pixels->Pixel2DClusterArray(p).at(imu);
        for ( auto const& pixel : cluster ) {
          cv::circle( cvimg, cv::Point(pixel.X(),pixel.Y()), 3, cv::Scalar(200,0,0), -1 );
        }
      }

      int nstopmu = ev_stopmu_pixels->Pixel2DClusterArray(p).size();
      for ( int imu=0; imu<nstopmu; imu++ ) {
        const larcv::Pixel2DCluster& cluster = ev_stopmu_pixels->Pixel2DClusterArray(p).at(imu);
        for ( auto const& pixel : cluster ) {
          cv::circle( cvimg, cv::Point(pixel.X(),pixel.Y()), 3, cv::Scalar(0,0,200), -1 );
        }
      }

      // draw proposed uncontained rois
      for ( auto const& roi : untagged_rois ) {
        larcv::draw_bb( cvimg, imgs_v.at(p).meta(), roi.BB().at(p), 0, 100, 100, 2 );        
      }

      // draw selected ROIs
      for ( auto const& roi : flash_matched_rois ){
        larcv::draw_bb( cvimg, imgs_v.at(p).meta(), roi.BB().at(p), 0, 200, 0, 2 );
      }

      std::stringstream ss;
      ss << "roiselected_n" << ientry << "_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
      cv::imwrite( ss.str(), cvimg );

    }

    // SAVE DATA TO OUTPUT FILE
    // ROIs, StopMu Tracks and clusters, ThruMu Tracks and Clusters, Truth ROI, Bad channels

    // save the image and bad channels
    larcv::EventImage2D* out_ev_images = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "tpc" );
    larcv::EventImage2D* out_ev_badchs = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "badchs" );
    for ( auto& img : imgs_v )
      out_ev_images->Append( img );
    for ( auto& img : badch_v )
      out_ev_badchs->Append( img );

    // save contained ROI
    larcv::EventROI*  out_ev_contained = (larcv::EventROI*)    dataco_output.get_larcv_data( larcv::kProductROI, "containedroi");
    larcv::EventPixel2D* out_ev_candidates = (larcv::EventPixel2D*) dataco_output.get_larcv_data( larcv::kProductPixel2D, "candidatepixels");
    out_ev_contained->Emplace( std::move(flash_matched_rois) );

    // truth quantities
    larcv::EventImage2D* out_ev_segmnt = (larcv::EventImage2D*)dataco_output.get_larcv_data( larcv::kProductImage2D, "segment");
    larcv::EventROI*     out_ev_roi    = (larcv::EventROI*)    dataco_output.get_larcv_data( larcv::kProductROI, "tpc");
    for ( auto& img : ev_seg->Image2DArray() )
      out_ev_segmnt->Append( img );
    out_ev_roi->Set( rois->ROIArray() );

    dataco_output.save_entry();

  }//end of event loop

  dataco_output.finalize();
  flash_matching.writeCalibTree();

  return 0;
}
