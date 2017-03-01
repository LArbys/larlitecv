#include "StopMuCluster.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {

  StopMuCluster::StopMuCluster( const StopMuClusterConfig& cfg ) : m_config(cfg) {
    // Constructor
    setVerbosity(m_config.verbosity);
  }

  void StopMuCluster::extractBaseClusters(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v, 
    const std::vector< std::vector< const larcv::Pixel2D* > >& endpts ) {

    // input checks
    if ( thrumu_v.size()!=img_v.size()) {
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of original and thrumu-tagged images are not the same.");
    }
    if ( img_v.size()!=m_config.pixel_thresholds.size() )
      throw std::runtime_error( "StopMuCluster::extractBaseClusters[Error] number of thresholds and images not the same.");

    // first mask out thrumu images
    m_masked_v.clear();
    m_img_v.clear();
    m_thrumu_v.clear();
    for ( size_t iimg=0; iimg<img_v.size(); iimg++ ) {
      const larcv::Image2D& img    = img_v.at(iimg);
      const larcv::Image2D& thrumu = thrumu_v.at(iimg);

      // we store pointers to these items for future reference
      m_img_v.push_back( &img );
      m_thrumu_v.push_back( &thrumu );

      // check that dimensions match
      if ( img.meta().rows()!=thrumu.meta().rows() || img.meta().cols()!=thrumu.meta().cols() ) {
        throw std::runtime_error("StopMuCluster::extractBaseClusters[Error] thrumu and orignal image dimensions do not match.");
      }

      larcv::Image2D masked( img.meta() );
      masked.paint(0.0);

      for (size_t r=0; r<img.meta().rows(); r++) {
        for (size_t c=0; c<img.meta().cols(); c++) {
          // skip if below threshold
          if ( img.pixel(r,c)<m_config.pixel_thresholds.at(iimg) ) 
            continue;
          // skip if tagged
          if ( thrumu.pixel(r,c)>0 )
            continue;
          masked.set_pixel(r,c,img.pixel(r,c));
        }
      }

      m_masked_v.emplace_back( std::move(masked) );


    }//end of image loop

    // unmask pixels around end points, so they can be included in the clusters
    for ( auto const& endpt : endpts ) {
      for (size_t p=0; p<m_masked_v.size(); p++ ) {
        larcv::Image2D& masked = m_masked_v.at(p);
        int row = (int)endpt.at(p)->Y();
        int col = (int)endpt.at(p)->X();
        for (int dr=-m_config.start_point_pixel_neighborhood; dr<m_config.start_point_pixel_neighborhood; dr++) {
          int r = row+dr;
          if ( r<0 || r>=masked.meta().rows()) continue;
          for (int dc=-m_config.start_point_pixel_neighborhood; dc<m_config.start_point_pixel_neighborhood; dc++) { 
            int c = col+dc;
            if ( c<0 || c>=masked.meta().cols() ) continue;
            masked.set_pixel(r,c,m_config.pixel_thresholds[p]+1);
          }
        }
      }
    }

    // use the masked image to form clusters
    m_untagged_clusters_v.clear();
    for (size_t p=0; p<m_masked_v.size(); p++) {
      untagged_cluster_info_t plane_cluster;
      plane_cluster.pixels = dbscan::extractPointsFromImage( m_masked_v.at(p), 0.5 );
      dbscan::DBSCANAlgo algo;
      plane_cluster.output = algo.scan( plane_cluster.pixels, m_config.dbscan_cluster_minpoints, m_config.dbscan_cluster_radius );
      for (size_t ic=0; ic<plane_cluster.output.clusters.size(); ic++) {
        dbscan::ClusterExtrema ex = dbscan::ClusterExtrema::FindClusterExtrema( ic, plane_cluster.output, plane_cluster.pixels );
        plane_cluster.extrema_v.emplace_back( std::move(ex) );
      }

      m_untagged_clusters_v.emplace_back( std::move(plane_cluster) );
    }

  }

  void StopMuCluster::saveClusterImageOCV( std::string filename ) {
    // for visual evaluation, we dump out various information used/constructed by this class
#ifdef USE_OPENCV
    std::vector<cv::Mat> cvimg_v = makeBaseClusterImageOCV();
    for (size_t p=0; p<cvimg_v.size(); p++) {
      std::stringstream ss;
      ss << filename << "_p" << p << ".png";
      cv::imwrite( ss.str(), cvimg_v.at(p) );
    }
#endif
  }

#ifdef USE_OPENCV
  std::vector<cv::Mat> StopMuCluster::makeBaseClusterImageOCV() {

    // we draw an image, that highlights the clusters, links between them, and interesting space points
    // we drawn an image per plane

    TRandom rand(1);
    std::vector<cv::Mat> cvimgs_v;

    for ( size_t p=0; p<m_img_v.size(); p++ ) {
      const larcv::Image2D& img = *(m_img_v.at(p));
      const untagged_cluster_info_t& cluster_info = m_untagged_clusters_v.at(p);

      // first make a CV image we can have fun with
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, m_config.pixel_thresholds[p], 100 );

      // color in the thrumu
      const larcv::Image2D& thrumu = *(m_thrumu_v.at(p));
      for (size_t r=0; r<thrumu.meta().rows(); r++) {
        for (size_t c=0; c<thrumu.meta().cols(); c++) {
          if ( thrumu.pixel(r,c)>0 ){
            cv::Vec3b& pix = cvimg.at<cv::Vec3b>( cv::Point(c,r) );
            pix[0] = 255;
            pix[1] = 0;
            pix[2] = 0;
          }
        }
      }

      // ok, we now label base clusters with colors!
      for (size_t icluster=0; icluster<cluster_info.output.clusters.size(); icluster++){
        if ( cluster_info.output.clusters.at(icluster).size()<m_config.dbscan_cluster_minpoints ) continue;
        // pick a color
        cv::Vec3b color;
        color[0] = (int)(rand.Uniform()*255);
        color[1] = (int)(rand.Uniform()*255);
        color[2] = (int)(rand.Uniform()*255);
        for (size_t ihit=0; ihit<cluster_info.output.clusters.at(icluster).size(); ihit++ ) {
          int hitidx = cluster_info.output.clusters.at(icluster).at(ihit);
          int x = cluster_info.pixels.at(hitidx)[0];
          int y = cluster_info.pixels.at(hitidx)[1];
          cvimg.at<cv::Vec3b>( cv::Point(x,y) ) = color;
        }
        // color in the extrema
        const dbscan::ClusterExtrema& ex = cluster_info.extrema_v.at(icluster);
        cv::circle(cvimg, cv::Point(ex.leftmost()[0],   ex.leftmost()[1]),   5, cv::Scalar(0,255,0),-1);
        cv::circle(cvimg, cv::Point(ex.topmost()[0],    ex.topmost()[1]),    5, cv::Scalar(255,255,0),-1);          
        cv::circle(cvimg, cv::Point(ex.rightmost()[0],  ex.rightmost()[1]),  5, cv::Scalar(0,255,255),-1);          
        cv::circle(cvimg, cv::Point(ex.bottommost()[0], ex.bottommost()[1]), 5, cv::Scalar(255,0,255),-1);
      }

      cvimgs_v.emplace_back( std::move(cvimg) );
    }

    return cvimgs_v;
  }
#endif  

}