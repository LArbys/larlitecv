#ifndef __ConfigBoundaryMuonTaggerAlgo__
#define __ConfigBoundaryMuonTaggerAlgo__

#include <vector>

// larcv
#include "Base/PSet.h"

namespace larlitecv {

  class ConfigBoundaryMuonTaggerAlgo {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can
      

  public:
    ConfigBoundaryMuonTaggerAlgo() {
      // defaults
      setdefaults();
    };
    ~ConfigBoundaryMuonTaggerAlgo() {};

    bool checkOK() { return true; }; // dummy for now
    
    std::vector<float> emptych_thresh; ///< pixel thresholds below which if wire stays, it is marked as empty (value per plane)
    std::vector<float> thresholds;    ///< pixel threshold to count as a hit
    std::vector<int>   neighborhoods; ///< columns before and after to check for hits
    std::vector<int>   edge_win_wires; ///
    std::vector<int>   edge_win_times;
    std::vector<float> edge_win_hitthresh;
    std::vector<int>   boundary_cluster_minpixels; ///< min number of pixels for dbscan clustering of boundary pixels
    std::vector<float>   boundary_cluster_radius;    ///< nearest neighbor radius for dbscan clustering of boundary pixels
    std::vector<float> astar_thresholds; //< passed to astar config
    std::vector<int>   astar_neighborhood; //< passed to astar config
    bool save_endpt_images;
    bool hitsearch_uses_badchs;
    int ticks_per_full_drift;
    std::vector<float> type_modifier;

    void setdefaults();
    
  };

  ConfigBoundaryMuonTaggerAlgo MakeConfigBoundaryMuonTaggerAlgoFromPSet( const larcv::PSet& pset );
  
}


#endif
