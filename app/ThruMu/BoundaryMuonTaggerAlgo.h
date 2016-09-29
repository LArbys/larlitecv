#ifndef __LLCV_BOUNDARY_MUON_TAGGER_ALGO__
#define __LLCV_BOUNDARY_MUON_TAGGER_ALGO__

#include <vector>
#include <string>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"

#include "BoundaryMatchArrays.h"

namespace larlitecv {

  class ConfigBoundaryMuonTaggerAlgo {
    // container for boundary muon tagger algo configuration parameters.
    // fill them anyway you can
      

  public:
    ConfigBoundaryMuonTaggerAlgo() {};
    ~ConfigBoundaryMuonTaggerAlgo() {};

    bool checkOK() { return true; }; // dummy for now

  public:
    
    std::vector<float> thresholds;    ///< pixel threshold to count as a hit
    std::vector<int>   neighborhoods; ///< columns before and after to check for hits
    
  };
    
  class BoundaryMuonTaggerAlgo {
    // the algo

  public:

    BoundaryMuonTaggerAlgo() {};
    BoundaryMuonTaggerAlgo( ConfigBoundaryMuonTaggerAlgo& config ) {
      _config = config; //copy
    };
    virtual ~BoundaryMuonTaggerAlgo() {};

    enum { kOK=0, kErr_NotConfigured, kErr_BadInput };

  public:
    
    void configure( ConfigBoundaryMuonTaggerAlgo& cfg ) { _config = cfg; };
    void run() {};
    int searchforboundarypixels( const std::vector< larcv::Image2D >& imgs, std::vector< larcv::Image2D >& boundarypixelimgs );

  protected:

    larlitecv::BoundaryMatchArrays m_matches;
    ConfigBoundaryMuonTaggerAlgo _config;

  };


};

#endif
