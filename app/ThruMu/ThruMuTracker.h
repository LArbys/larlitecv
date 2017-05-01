#ifndef __THRUMU_TRACKER__
#define __THRUMU_TRACKER__

// stdlib
#include <vector>

// larcv
#include "DataFormat/Image2D.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "ThruMuTrackerConfig.h"
#include "Linear3DChargeTagger.h"

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

    struct AStarTaggerInfo {
      bool isgood;
      bool goal_reached;
      int nbad_nodes;
      int total_nodes;
      float frac_bad;
      AStarTaggerInfo() {
        isgood = false;
	goal_reached = false;
	nbad_nodes = -1;
	total_nodes = -1;
	frac_bad = 0.;
      };
    };


    void runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
		  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
		  std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters );

    larlitecv::BMTrackCluster3D runLinearChargeTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                    const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                    LinearTaggerInfo& result_info );

    larlitecv::BMTrackCluster3D runAStarTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                   const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                   AStarTaggerInfo& result_info );


    // for astar, compressed image containers
    std::vector< larcv::Image2D > m_img_compressed_v;
    std::vector< larcv::Image2D > m_badch_compressed_v;
    
    
  };

}

#endif
