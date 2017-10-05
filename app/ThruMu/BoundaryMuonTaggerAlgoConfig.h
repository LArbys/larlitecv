#ifndef __BoundaryMuonTaggerAlgoConfig__
#define __BoundaryMuonTaggerAlgoConfig__

#include <vector>

// larcv
#include "Base/PSet.h"

namespace larlitecv {

  class BoundaryMuonTaggerAlgoConfig {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can

  public:
    BoundaryMuonTaggerAlgoConfig() {
      // defaults
      setdefaults();
    };
    ~BoundaryMuonTaggerAlgoConfig() {
    };

    bool checkOK() { return true; }; // dummy for now
    void print();

    std::vector<float> thresholds;    ///< pixel threshold to count as a hit
    std::vector<int>   neighborhoods; ///< columns before and after to check for hits
    std::vector<int>   boundary_cluster_minpixels; ///< min number of pixels for dbscan clustering of boundary pixels
    std::vector<float> boundary_cluster_radius;    ///< nearest neighbor radius for dbscan clustering of boundary pixels
    bool save_endpt_images;
    bool hitsearch_uses_badchs;
    std::vector<float> type_modifier;
    int fKernelRadius;
    int verbosity;

    void setdefaults();

  };

  // larcv/larlitecv interface
  BoundaryMuonTaggerAlgoConfig MakeBoundaryMuonTaggerAlgoConfigFromPSet( const larcv::PSet& pset );

}


#endif
