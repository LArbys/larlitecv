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

// larlitecv
#include "ContainedROIConfig.h"

namespace larlitecv {


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
    	float total_charge;
    };

    // Primary Method
    std::vector<larcv::ROI> SelectROIs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v,
    	const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v, 
    	const std::vector<larlite::event_opflash*>& opflashes_v );


    // Supporting Methods
		larcv::Image2D CollectedUntaggedPixels( const larcv::Image2D& img, const larcv::Image2D& thrumu, const larcv::Image2D& stopmu);
		std::vector< analyzed_cluster_t > AnalyzeClusters( const untagged_cluster_info_t& clusters_info, const larcv::Image2D& img );
		float CalculateMatchLikelihood( const std::vector<analyzed_cluster_t>& clusters, const std::vector<larcv::Image2D>& img_v, std::vector<float>& ll_components );
    void MatchClusters( const std::vector<analyzed_cluster_t>& uplane, const std::vector<analyzed_cluster_t>& vplane, const std::vector<analyzed_cluster_t>& yplane, 
			const std::vector<larcv::Image2D>& img_v,  std::vector< std::vector<analyzed_cluster_t> >& cluster_combos, 
			std::map< int, float >& combo_likelihoods, std::map<int,float>& combo_best_triarea);

    int m_run;
    int m_subrun;
    int m_event;
    int m_entry;
  };
 

}



#endif
