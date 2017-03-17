#ifndef __ConfigBoundaryMuonTaggerAlgo__
#define __ConfigBoundaryMuonTaggerAlgo__

#include <vector>

// larcv
#include "Base/PSet.h"
#include "AStar3DAlgo.h"
#include "Linear3DChargeTagger.h"

namespace larlitecv {

  class ConfigBoundaryMuonTaggerAlgo {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can
      

  public:
    ConfigBoundaryMuonTaggerAlgo() {
      // defaults
      setdefaults();
    };
    ~ConfigBoundaryMuonTaggerAlgo() {
    };

    bool checkOK() { return true; }; // dummy for now
    void print();
    
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
    std::vector<int>   tag_neighborhood; //< neighborhood around track path that will be tagged as thrumu
    bool save_endpt_images;
    bool hitsearch_uses_badchs;
    int ticks_per_full_drift;
    std::vector<float> type_modifier;
    int verbosity;
    AStar3DAlgoConfig    astar_cfg;
    Linear3DChargeTaggerConfig linear3d_cfg;
    int linear3d_min_tracksize;
    float linear3d_min_goodfraction;
    float linear3d_min_majoritychargefraction;
    float astar3d_min_goodfrac;
    float astar3d_min_majfrac;

    void setdefaults();
    const AStar3DAlgoConfig& getAStarConfig() { return astar_cfg; };
    const Linear3DChargeTaggerConfig& getLinear3DConfig() { return linear3d_cfg; };
    
  };

  ConfigBoundaryMuonTaggerAlgo MakeConfigBoundaryMuonTaggerAlgoFromPSet( const larcv::PSet& pset );
  
}


#endif
