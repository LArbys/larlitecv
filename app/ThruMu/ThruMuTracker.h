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
#include <algorithm>
#include <set>

// larcv
#include "DataFormat/Image2D.h"

// larlite
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "OpT0Finder/Algorithms/TimeCompatMatch.h"
#include "OpT0Finder/Base/BaseAlgorithm.h"
#include "OpT0Finder/Base/BaseProhibitAlgo.h"
#include "FhiclLite/FhiclLiteUtilFunc.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "ThruMuTrackerConfig.h"
#include "Linear3DChargeTagger.h"
#include "FlashMatchInterface/GeneralFlashMatchAlgoConfig.h"
#include "FlashMatchInterface/GeneralFlashMatchAlgo.h"

namespace larlitecv {

  class ThruMuTracker {

  public:

    ThruMuTracker( const ThruMuTrackerConfig& config );
    virtual ~ThruMuTracker() {};

    // Add two new arguments onto the start of the function:
    // fcllite::PSet - This is the PSet necessary for configuring the 'FlashMatchManager' in the 'flashMatchTracks' function.
    // 'flash_match_config' - this is the set of configuration parameters for the 'GeneralFlashMatchAlgo' class, of type 'GeneralFlashMatchAlgoConfig'.

    // Add two new arguments onto the end of this function:
    // flash_idx_v: This is the vector of the indices that correspond to the endpoints contained in the 'spacepts' argument.
    //  Note that this vector
    //  (1) assumes that the two flash producers are 'simpleFlashBeam' and 'simpleFlashCosmic' and
    //  (2) fills the flash indices in a running, one-dimensional list starting with the flashes produced by 'simpleFlashBeam'
    //      before continuing onto the flashes produced by 'simpleFlashCosmic'.
    void makeTrackClusters3D( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector<larcv::Image2D>& img_v,
			      const std::vector<larcv::Image2D>& badchimg_v, const std::vector< const BoundarySpacePoint* >& spacepts,
			      std::vector< larlitecv::BMTrackCluster3D >& trackclusters,
			      std::vector< larcv::Image2D >& tagged_v, std::vector<int>& used_endpoints_indices, const std::vector< larlite::event_opflash* >& opflashsets);

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

    // Declare the structure for the information concerning the candidate tracks for a flash.
    struct FlashMatchCandidate {
      int idx; // the track's original index in the qcluster
      float chi2;  // the score between the track and the flash
    };

    // Declare a structure for the information in finding a unique track/flash match.
    // This is identical to 'FlashMatchCandidate' except that it contains 'spot_in_ranking', which is where this qcluster stands in the qcluster ranking for the flash.
    // These pieces of information will be contained in a vector, with the spot in the vector corresponding to the location of this flash in the 'all_opflash_candidate_list' vector.
    struct UniqueTrackFlashMatch {
      int idx;
      float chi2;
      int spot_in_ranking;
    };

    void runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg, const std::vector< const BoundarySpacePoint* >& spacepts,
		  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v, std::vector<larcv::Image2D>& tagged_v,
		  std::vector<int>& used_endpoints_indices, std::vector<larlitecv::BMTrackCluster3D>& trackclusters, std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v,
		  std::vector< std::vector< larlitecv::BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v, std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v );

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
    void flashMatchTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< const BoundarySpacePoint* >& spacepts, const std::vector< larlite::event_opflash* >& opflash_v, 
					  std::vector< BMTrackCluster3D >& trackclusters, std::vector< std::vector< BoundarySpacePoint >  >& impossible_match_endpoint_v, 
			                  std::vector< BoundaryFlashIndex >& already_matched_flash_v, std::vector< int >& well_matched_tracks_idx_v, const int& num_of_tracks_added_in_pass, 
			                  std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v, 
			                  std::vector< std::vector< larlitecv::BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v, std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v, 
			                  bool anode_and_cathode_only );
		
    void flashMatchAC( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< flashana::QCluster_t >& qcluster_vector,
				      const std::vector< BoundaryFlashIndex >& boundary_flash_index_vector, const std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v_from_pass,
				      const std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v_from_pass,
				      std::vector< std::vector< BoundarySpacePoint > >& impossible_match_endpoint_v, const std::vector< std::vector< BoundarySpacePoint > > & track_endpoint_v_from_pass,
		                      std::vector < int >& well_matched_tracks_idx_v, std::vector< BoundaryFlashIndex >&  already_matched_flash_v );

    void flashMatchYZFaceTracks( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< larlitecv::BMTrackCluster3D >& trackclusters, 
				 const std::vector < flashana::QCluster_t >& qcluster_vector, const std::vector< BoundaryFlashIndex >& boundary_flash_index_vector, 
				 std::vector< std::vector< BoundarySpacePoint > >& impossible_match_endpoint_v, const std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v_from_pass, 
				 std::vector< int >& well_matched_tracks_idx_v, std::vector< BoundaryFlashIndex >& already_matched_flash_v, bool entire_event );

    void rankTrackFlashMatches( GeneralFlashMatchAlgoConfig& flash_match_config, const std::vector< flashana::QCluster_t >& qclusters_being_checked, const std::vector< int >& well_matched_tracks_idx_v, 
				larlite::opflash opflash, int opflash_producer, std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list );
    
    void orderInAscendingChi2Order( const std::vector< ThruMuTracker::FlashMatchCandidate >& opflash_track_match_list, 
						    std::vector< ThruMuTracker::FlashMatchCandidate >& ordered_opflash_track_match_list );

    void findUniqueTrackFlashMatch( const std::vector< std::vector< ThruMuTracker::FlashMatchCandidate > >& all_opflash_candidate_list, std::vector< float >& best_chi2_v, 
		       std::vector< int >& corresponding_qcluster_idx_v );

    bool isSameQClusterMatchedToDifferentFlashes( const std::vector< ThruMuTracker::UniqueTrackFlashMatch >& unique_track_flash_match_v );
      
    void findQClusterFlashMatchWithWorseChi2( const std::vector< ThruMuTracker::UniqueTrackFlashMatch >& unique_track_flash_match_v, int& unique_track_flash_match_idx, 
					      ThruMuTracker::UniqueTrackFlashMatch& duplicate_flash_match_with_worse_chi2 );

    void sortOutBadTracks( std::vector < larlitecv::BMTrackCluster3D >& trackclusters, const std::vector< int >& well_matched_tracks_idx_v, std::vector< int >& tracks_per_pass, int tracks_under_consideration, bool single_pass );
		


    // for astar, compressed image containers
    std::vector< larcv::Image2D > m_img_compressed_v;
    std::vector< larcv::Image2D > m_badch_compressed_v;


  };

}

#endif
