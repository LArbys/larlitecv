#include <iostream>
#include <string>
#include <vector>

#include "DataFormat/EventImage2D.h"
#include "DataFormat/Image2D.h"

#include "Base/DataCoordinator.h"

#include <opencv2/opencv.hpp>
#include "CVUtil/CVUtil.h"

#include "ClusterGroupAlgo.h"

int main( int nargs, char** argv ) {

  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( "test_larcv.root", "larcv" );
  dataco.configure("croi3d.cfg", "StorageManager", "IOManager", "CROI3D" );
  dataco.initialize();

  dataco.goto_entry(0,"larcv");

  larcv::EventImage2D* ev_img    = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimgs" );
  larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "marked3d" );
  larcv::EventImage2D* ev_gapchs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "gapchs" );  
  larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "stopmu" );

  const std::vector<larcv::Image2D>& img_v    = ev_img->Image2DArray();
  const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
  const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();
  const std::vector<larcv::Image2D>& gapchs_v = ev_gapchs->Image2DArray();

  // first attempt at image subtraction. what does that look like?
  std::vector< larcv::Image2D > subimg_v;
  for ( size_t p=0; p<img_v.size(); p++) {
    larcv::Image2D sub( img_v.at(p) );
    for ( size_t r=0; r<sub.meta().rows(); r++ ) {
      for ( size_t c=0; c<sub.meta().cols(); c++ ) {
      	float val = sub.pixel(r,c);
      	if ( thrumu_v.at(p).pixel(r,c)>0 || stopmu_v.at(p).pixel(r,c)>0 ) {
      	  val = 0.0;
      	}
      	sub.set_pixel(r,c,val);
      }
    }
    cv::Mat cvimg = larcv::as_mat_greyscale2bgr( sub, 10, 60 );
    std::stringstream ss;
    ss << "subimg_p" << p << ".jpg";
    cv::imwrite( ss.str(), cvimg );
    subimg_v.emplace_back( std::move(sub) );
  }

  larlitecv::ClusterGroupAlgoConfig config;
  larlitecv::ClusterGroupAlgo cluster_algo(config);

  cluster_algo.MakeClusterGroups( img_v, gapchs_v, subimg_v );

  return 0;
}
