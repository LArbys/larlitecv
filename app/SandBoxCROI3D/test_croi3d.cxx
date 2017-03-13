#include <iostream>
#include <string>
#include <vector>

#include "DataFormat/EventImage2D.h"
#include "DataFormat/Image2D.h"

#include "Base/DataCoordinator.h"

#include <opencv2/opencv.hpp>
#include "CVUtil/CVUtil.h"

#include "ClusterGroupAlgo.h"
#include "ClusterGroupMatchingAlgo.h"

#include "TRandom.h"

int main( int nargs, char** argv ) {

  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( "test_larcv.root", "larcv" );
  dataco.configure("croi3d.cfg", "StorageManager", "IOManager", "CROI3D" );
  dataco.initialize();

  dataco.goto_entry(1,"larcv");

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
      	if ( val<10.0 || thrumu_v.at(p).pixel(r,c)>0 || stopmu_v.at(p).pixel(r,c)>0 ) {
      	  val = 0.0;
      	}
      	sub.set_pixel(r,c,val);
      }
    }

    //cv::Mat cvimg = larcv::as_mat_greyscale2bgr( sub, 10, 60 );
    //std::stringstream ss;
    //ss << "subimg_p" << p << ".jpg";
    //cv::imwrite( ss.str(), cvimg );

    subimg_v.emplace_back( std::move(sub) );
  }

  // form tagged image
  std::vector< larcv::Image2D > tagged_v;
  for ( size_t p=0; p<img_v.size(); p++) {
    larcv::Image2D tagged( img_v.at(p).meta() );
    tagged.paint(0.0);
    for ( size_t r=0; r<tagged.meta().rows(); r++ ) {
      for ( size_t c=0; c<tagged.meta().cols(); c++ ) {
      	if ( thrumu_v.at(p).pixel(r,c)>0 || stopmu_v.at(p).pixel(r,c)>0 )
          tagged.set_pixel(r,c,255);
      }
    }
    tagged_v.emplace_back( std::move(tagged) );
  }  

  larlitecv::ClusterGroupAlgoConfig config;
  config.dbscan_cluster_minpoints = 20;
  config.alldir_max_link_dist = 30.0;
  config.max_link_distance = 300.0;
  config.min_link_cosine = 0.9;
  config.single_cluster_group_min_npoints = 100;
  larlitecv::ClusterGroupAlgo cluster_algo(config);

  std::vector< larlitecv::PlaneClusterGroups > plane_groups = cluster_algo.MakeClusterGroups( img_v, gapchs_v, tagged_v );
  std::vector< cv::Mat > cvimg_cluster;


  larlitecv::ClusterGroupMatchingAlgo matcher;
  std::vector<larlitecv::ChargeVolume> vols = matcher.MatchClusterGroups( subimg_v, plane_groups );

  // image the output of clustergroup matcher
  std::vector<cv::Mat> cv_matched;
  for ( size_t p=0; p<img_v.size(); p++) {
    cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img_v.at(p), 10, 60 );
    cv_matched.emplace_back( std::move(cvimg) );
  }
  matcher.labelCVImageWithMatchedClusters( cv_matched, img_v, vols, 0.75 );
  for ( size_t p=0; p<cv_matched.size(); p++) {
    cv::Mat& cvimg = cv_matched.at(p);
    std::stringstream ss;
    ss << "clustermatch_p" << p << ".jpg";
    imwrite( ss.str(), cvimg );
  }

  return 0;
}
