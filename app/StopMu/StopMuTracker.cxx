#include "StopMuTracker.h"

#include <stdio.h>
#include <time.h>
#include <sstream>

#ifndef __CINT__
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif

#include "StopMuSkeleton.h"

namespace larlitecv {

  StopMuTracker::StopMuTracker( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v ) {
				
    
    // skeletonize image
    StopMuSkeleton skeleton_op;

    for (int p=0; p<3; p++) {
      larcv::Image2D skel = skeleton_op.skeletonize( img_v.at(p), 10.0, 3 );
      skel_v.emplace_back( std::move(skel) );
    }
    
    // cluster skeleton pixel, but mask thru-mu pixels
    time_t cluster_start = time(NULL);
    
    for (int p=0; p<3; p++) {
      dbscan::dbPoints data;
      const larcv::Image2D& skelimg = skel_v.at(p);
      for (int r=0; r<skelimg.meta().rows(); r++) {
	for (int c=0; c<skelimg.meta().cols(); c++) {
	  if ( skelimg.pixel(r,c)==0 ) continue; // not skeleton
	  if ( thrumu_v.at(p).pixel(r,c)>0 ) continue; // mask thrumu
	  std::vector<double> point(2);
	  point[0] = c; // X 
	  point[1] = r; // Y
	  data.emplace_back( point );
	}
      }
      dbscan::DBSCANAlgo dbalgo;
      dbscan::dbscanOutput cluster = dbalgo.scan( data, 5, 5.0 );
      imghits.emplace_back( data );
      clusters.emplace_back( cluster );
      std::cout << "number of clusters on plane " << p << ": " << clusters.at(p).clusters.size() << std::endl;
    }//loop over clusters
    
    time_t cluster_finished = time(NULL);
    //double dt_clusters = difftime(cluster_start,cluster_finished);
    double dt_clusters = cluster_finished-cluster_start;
    std::cout << "clustered in " << dt_clusters << " seconds." << std::endl;

  }

  void StopMuTracker::trackStopMu( const std::vector< std::vector<int> >& start2d, const std::vector< std::vector<float> >& start_dir2d,
				   const std::vector< float >& start_pos3d, const std::vector<float>& start_dir3d ) {
    // find matching cluster
    std::vector<int> clusterid;

    const larcv::ImageMeta& meta = skel_v.at(0).meta();
    larcv::Image2D img_cluster( meta );
    img_cluster.paint(0.0);

    for (int p=0; p<3; p++) {
      std::vector<double> testpoint(2);
      testpoint[0] = start2d.at(p)[0]; // X
      testpoint[1] = start2d.at(p)[1]; // Y
      std::cout << "plane " << p << " testpoint: (" << testpoint[0] << "," << testpoint[1] << ")" << std::endl;
      int match = clusters.at(p).findMatchingCluster( testpoint, imghits.at(p), 5.0 );
      clusterid.push_back(match);
      std::cout << "plane=" << p << " matching cluster index=" << match << std::endl;

      // we make a sorted list of pixels by distance
      // we also need an initial direction
      hitlists[p] = new Hit2DList;

      for (int ihit=0; ihit<clusters.at(p).clusters.at(match).size(); ihit++) {
	int hitidx = clusters.at(p).clusters.at(match).at(ihit);
	int x = imghits.at(p).at(hitidx)[0];
	int y = imghits.at(p).at(hitidx)[1];

	// dir from start to point
	std::vector<float> dir(2);
	dir[0] = x-start2d.at(p)[0];
	dir[1] = y-start2d.at(p)[1];
	float norm = sqrt( dir[0]*dir[0] + dir[1]*dir[1] );
	for (int i=0; i<2; i++) dir[i] /= norm;
	
	float cosine = 0.;
	for (int i=0; i<2; i++) cosine += dir[i]*start_dir2d.at(p)[i];

	Hit2D hit;
	hit[0] = x;
	hit[1] = y;
	if ( cosine>0 )
	  hit.distance = norm;
	else
	  hit.distance = -norm;
	if ( cosine>0 ) {
	  hitlists[p]->emplace(std::move(hit));
	  hitlists[p]->sort();
	}
	
	img_cluster.set_pixel((int)y,(int)x,250.0);
      }

      
      std::cout << "plane " << p << " number of hits=" << hitlists[p]->size() << std::endl;
      for (int ihit=0; ihit<(int)hitlists[p]->size(); ihit++) {
	std::cout << " [#" << ihit << "] (" << hitlists[p]->at(ihit)[0] << "," << hitlists[p]->at(ihit)[1] << ")"
		  << " " << hitlists[p]->at(ihit).distance << std::endl; 
      }
      
    }//end of loop over planes
    
    cv::Mat imgmat = larcv::as_mat( img_cluster );
    std::stringstream ss;
    ss << "baka.jpg";
    cv::imwrite( ss.str().c_str(), imgmat );
    
  }

}
