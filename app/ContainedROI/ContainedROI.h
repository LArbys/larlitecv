#ifndef __CONTAINED_ROI__
#define __CONTAINED_ROI__

/**
   This class is responsible for choosing candidate contained clusters.
**/

#include <vector>

// larcv
#include "DataFormat/Image2D.h"
#include "DataFormat/ROI.h"
#include "DataFormat/Pixel2D.h"
#include "dbscan/DBSCANAlgo.h"

// larlite
#include "DataFormat/opflash.h"

namespace larlitecv {

	class ContainedROIConfig {
	public:
		ContainedROIConfig() {};
		virtual ~ContainedROIConfig() {};

		float pixel_threshold;
		int min_cluster_size;

	};

  class ContainedROI {
  public:

    ContainedROI( const ContainedROIConfig& config );
    virtual ~ContainedROI() {};

    ContainedROIConfig m_config;

    struct untagged_cluster_info_t {
    	::dbscan::dbPoints pixels;
    	::dbscan::dbscanOutput output;
    };

    struct analyzed_cluster_t {
    	int cluster_idx;
    	std::vector< larcv::Pixel2D > extrema_pts;
    	std::vector<float> mean;
    };

    // Primary Method
    std::vector<larcv::ROI> SelectROIs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v,
    	const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v, 
    	const larlite::event_opflash& opflash_v );


    // Supporting Methods
		larcv::Image2D CollectedUntaggedPixels( const larcv::Image2D& img, const larcv::Image2D& thrumu, const larcv::Image2D& stopmu);
		std::vector< analyzed_cluster_t > AnalyzeClusters( const untagged_cluster_info_t& clusters_info );
    
  };


}



#endif
