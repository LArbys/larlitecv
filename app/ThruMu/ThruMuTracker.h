#ifndef __THRUMU_TRACKER__ 
#define __THRUMU_TRACKER__ 

// stdlib
#include <vector>

// larcv
#include "DataFormat/Image2D.h"

// larlitecv
#include "BoundarySpacePoint.h"
#include "BMTrackCluster3D.h"
#include "ThruMuTrackerConfig.h"

namespace larlitecv {

  class ThruMuTracker {

  public:

    ThruMuTracker( const ThruMuTrackerConfig& config );
    virtual ~ThruMuTracker() {};


    void makeTrackClusters3D( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
			      const std::vector< const BoundarySpacePoint* >& spacepts,
			      std::vector< larlitecv::BMTrackCluster3D >& trackclusters, 
			      std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices);    

  protected:

    ThruMuTrackerConfig m_config;

    // we use this class to communicate how good the linear track hypothesis was
    struct LinearTaggerInfo {
      bool isgood;
      int numpts;
      float goodfrac;
      float majfrac;
      LinearTaggerInfo() { 
	isgood=false;
	numpts = 0;
	goodfrac = 0.;
	majfrac = 0.;
      };
    };
    

    void runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
		  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
		  std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters );

    larlitecv::BMTrackCluster3D runLinearChargeTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg, 
						       const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
						       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
						       LinearTaggerInfo& result_info );




  };

}

#endif
