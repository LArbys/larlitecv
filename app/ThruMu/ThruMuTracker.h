#ifndef __THRUMU_TRACKER__
#define __THRUMU_TRACKER__

/* ====================================================
 *  ThruMuTracker
 *
 *  We use a number of different tracking
 *  algorithms to solve the thrumu-problem.
 *
 *
 * ==================================================== */

// stdlib
#include <vector>

// larcv
#include "DataFormat/Image2D.h"

// larlite
#include "DataFormat/opflash.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "ThruMuTrackerConfig.h"
#include "Linear3DChargeTagger.h"
#include "GeneralFlashMatchAlgo/GeneralFlashMatchAlgoConfig.h"
#include "GeneralFlashMatchAlgo/GeneralFlashMatchAlgo.h"

namespace larlitecv {

  class ThruMuTracker {

  public:

    ThruMuTracker( const ThruMuTrackerConfig& config );
    virtual ~ThruMuTracker() {};

    // Add an argument onto the start of the function:
    // 'flash_config' - this is the set of configuration parameters for the 'GeneralFlashMatchAlgo' class, of type 'GeneralFlashMatchAlgoConfig'.

    // Add two new arguments onto the end of this function:
    // flash_idx_v: This is the vector of the indices that correspond to the endpoints contained in the 'spacepts' argument.  Note that this vector (1) assumes that the two flash producers are 'simpleFlashBeam' and 'simpleFlashCosmic' and (2) fills the flash indices in a running, one-dimensional list starting with the flashes produced by 'simpleFlashBeam' before continuing onto the flashes produced by 'simpleFlashCosmic'.
    void makeTrackClusters3D( GeneralFlashMatchAlgoConfig& flash_config, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
					     const std::vector< const BoundarySpacePoint* >& spacepts,
					     std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
					     std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices, const std::vector< larlite::event_opflash* >& opflashsets,
					     std::vector< int > & flash_idx_v, std::vector< int >& flash_producer_idx_v );

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

    struct FoxTrotExtenderInfo {
      bool isgood;
      FoxTrotExtenderInfo() {
        isgood = false;  //< extended back
      };
    };


    void runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
		  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
		  std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters, const std::vector< int >& flash_idx_v,
                               const std::vector< int >& flash_producer_idx_v, std::vector< int > track_endpoint_flash_idx_v, std::vector< int > track_endpoint_flash_producer_idx_v );

    larlitecv::BMTrackCluster3D runLinearChargeTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                    const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                    LinearTaggerInfo& result_info );

    larlitecv::BMTrackCluster3D runAStarTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                   const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                   AStarTaggerInfo& result_info );

    void runFoxTrotExtender( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg, std::vector<std::vector<double> >& track,
                  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
                  ThruMuTracker::FoxTrotExtenderInfo& result_info );

    // Add the new function for flash-matching the tracks in the ThruMu stage of reconstruction.
    void flashMatchTracks( GeneralFlashMatchAlgoConfig& flash_config, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& tagged_v, 
			   const std::vector< const BoundarySpacePoint* >& spacepts, const std::vector< larlite::event_opflash* >& opflash_v, std::vector< BMTrackCluster3D >& trackclusters, 
			   std::vector< int >& impossible_match_endpoints, std::vector< int >& already_matched_flash_idx, const int& num_of_tracks_added_in_pass, 
			   std::vector< int >& track_endpoint_flash_idx_v, std::vector< int >& track_endpoint_flash_producer_idx_v );


    // for astar, compressed image containers
    std::vector< larcv::Image2D > m_img_compressed_v;
    std::vector< larcv::Image2D > m_badch_compressed_v;


  };

}

#endif
