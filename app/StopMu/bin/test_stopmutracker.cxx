#include <iostream>
#include <sstream>

// config/storage
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

#ifndef __CINT__
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "CVUtil/CVUtil.h"

// larlitecv
#include "ThruMu/EmptyChannelAlgo.h"

// larlite/stopmu
#include "StopMuAlgoTypes.h"
#include "StopMuStart.h"
#include "StopMuTracker.h"
#include "StopMuSkeleton.h"


int main( int nargs, char** argv ) {

  std::cout << "Test the stop mu tracker." << std::endl;

  // config file
  std::string cfg_file = "kfstop.cfg";

  // data coordinator: load an example data file
  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( "output_larcv_testextbnb_001.root", "larcv" );
  dataco.configure( cfg_file, "StorageManager", "IOManager", "KFStopMu" );
  dataco.initialize();

  // go to some entry
  dataco.goto_entry(0, "larcv");

  // STOP MU START PARAMETERS
  float fThreshold = 10.0;
  int rneighbor = 10;
  int cneighbor = 10;

  // for the test, we target a top-passing, stop muon
  // at time around tick 3880
  larcv::EventImage2D* imgs            = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
  larcv::EventImage2D* marked3d        = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
  larcv::EventPixel2D* top_spacepoints = (larcv::EventPixel2D*)dataco.get_larcv_data( larcv::kProductPixel2D, "topspacepts" );
  
  // make the bad channel image
  larlitecv::EmptyChannelAlgo emptyalgo;
  larcv::EventChStatus* ev_status      = (larcv::EventChStatus*)dataco.get_larcv_data( larcv::kProductChStatus, "tpc" );
  std::vector< larcv::Image2D > badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status ); //< fill this out!
  std::cout << "number of bad ch imgs: " << badch_v.size() << std::endl;
  
  const std::vector<larcv::Image2D>& img_v    = imgs->Image2DArray();
  const std::vector<larcv::Image2D>& thrumu_v = marked3d->Image2DArray();
  const larcv::ImageMeta& meta = img_v.at(0).meta();
  
  int nendpts = top_spacepoints->Pixel2DArray(0).size();
  std::vector<larcv::Pixel2D> start;
  for (int i=0; i<nendpts; i++) {
    const larcv::Pixel2D& pix = top_spacepoints->Pixel2DArray(0).at(i);
    float tick = meta.pos_y( pix.Y() );
    std::cout << "top endpoint #" << i << ": tick=" << tick << std::endl;
    //if ( tick>3850 && tick<3950 ) { // test point A
    if ( tick>6900 && tick<6910 ) { // test point B
      std::cout << "Found test start point:tick= " << tick << std::endl;
      for (int p=0; p<3; p++) {
	larcv::Pixel2D copy( top_spacepoints->Pixel2DArray(p).at(i) );
	start.emplace_back( copy );
      }
    }
  }
  
  std::cout << "start point: " << start.size() << std::endl;

  // --------------------------------------------
  // starting point algorithm
  larlitecv::StopMuStart start_finder_algo;
  start_finder_algo.setVerbose(1);
  std::vector<float> start_spacepoint;
  std::vector< std::vector<float> > start_dir2d;
  std::vector<float> start_dir3d;
  start_finder_algo.getStartDirectionV( img_v, badch_v, start, rneighbor, cneighbor, fThreshold, start_spacepoint, start_dir2d, start_dir3d );

  // --------------------------------------------
  // pass start info and image into stopmu tracker
  // translate input
  std::vector< std::vector<int> > start2d_pos;
  for (int p=0; p<3; p++) {
    std::vector<int> pos(2);
    pos[0] = start.at(p).X();
    pos[1] = start.at(p).Y();
    start2d_pos.emplace_back( pos );
  }
  larlitecv::StopMuTracker sttracker( img_v, thrumu_v );
  larlitecv::Step3D start_track;
  sttracker.trackStopMu( start2d_pos, start_dir2d, start_spacepoint, start_dir3d, start_track );

  // let's mark up an image
  std::vector<larcv::Image2D> out_imgs;
  for (int p=0; p<3; p++) {
    larcv::Image2D img( img_v.at(p) );
    for (int r=0; r<img.meta().rows(); r++) {
      for (int c=0; c<img.meta().cols(); c++) {
	if ( img.pixel(r,c)<10 )
	  img.set_pixel(r,c,0);
	else if ( img.pixel(r,c)>200 )
	  img.set_pixel(r,c,200);
      }
    }
    out_imgs.emplace_back( std::move(img) );
  }
  larlitecv::Step3D& current_step = start_track;
  int nsteps = 0;
  do {
    if ( current_step.planepositions.size()!=3 ) {
      nsteps++;
      if ( !current_step.isEnd() )
	current_step = current_step.GetNext();
      continue;
    }
    //std::cout << "saving step: " << nsteps << std::endl;
    for (int p=0; p<3; p++) {
      out_imgs.at(p).set_pixel( current_step.planepositions.at(p)[1], current_step.planepositions.at(p)[0], 255 );
    }
    if ( !current_step.isEnd() )
      current_step = current_step.GetNext();
    nsteps++;
  } while ( !current_step.isEnd() );

  for (int p=0; p<3; p++) {
    const larcv::Image2D& outimg = out_imgs.at(p);
    cv::Mat imgmat = larcv::as_mat_1FC( outimg );
    std::stringstream ss;
    ss << "stmu_tracker_out_plane" << p << ".jpg";
    cv::imwrite( ss.str(), imgmat );
  }
  
  return 0;
}
