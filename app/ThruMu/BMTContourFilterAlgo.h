#ifndef __BMT_CONTOUR_FILTER_H__
#define __BMT_CONTOUR_FILTER_H__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

#include "TaggerContourTools/ContourCluster.h"
#include "TaggerContourTools/ContourShapeMeta.h"

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif


namespace larlitecv {

  
  class BMTContourFilterAlgo {
  public:

    BMTContourFilterAlgo();
    virtual ~BMTContourFilterAlgo();


    bool buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
		       const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
		       const float max_dist2contour );

    bool extendClusterGroup( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
    			     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v );
    
    /* std::vector<float> calculateContourIntersection( const std::vector< cv::Point >& cnt1, const std::vector< cv::Point >& cnt2 ); */

    /* bool buildContourGraph( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v, */
    /* 			    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ); */

    int getIndexOfContainingContour( const int row, const int col, const std::vector<ContourShapeMeta>& contours_v, int min_cluster_size, float dist_tolerance );
    
    bool ratchetCluster( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
    			 const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v );

    bool isPointInContour( const std::vector<cv::Point>& imgpt, const std::vector<larcv::Image2D>& img_v,
			   const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
			   const float max_dist2contour,
			   ContourCluster& outcluster );
    
    bool analyzeSeedContours( const std::vector<cv::Point>& imgpt, const std::vector<larcv::Image2D>& img_v,
			      const float max_dist2edge,
			      const ContourCluster& seedcluster );

    // accessable intermediate variables for debugging
    int num_valid_pos;
    int num_combos;
    std::vector< int > min_otherplane_v;
    std::vector< int > max_otherplane_v; 
    std::vector< std::vector<float> > max_poszy_v;
    std::vector< std::vector<float> > min_poszy_v;    
    std::vector<bool> cuts;
    
  };

}


#endif
