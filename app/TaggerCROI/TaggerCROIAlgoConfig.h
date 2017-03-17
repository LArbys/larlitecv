#ifndef __TAGGER_CROI_ALGO_CONFIG_H__
#define __TAGGER_CROI_ALGO_CONFIG_H__

// larcv
#include "Base/PSet.h"

// larlitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"

namespace larlitecv {

  class TaggerCROIAlgoConfig {
  public:
    TaggerCROIAlgoConfig() {};
    virtual ~TaggerCROIAlgoConfig() {};

    larlitecv::ConfigBoundaryMuonTaggerAlgo sidetagger_cfg;
		larlitecv::FlashMuonTaggerConfig        flashtagger_cfg;

    static TaggerCROIAlgoConfig makeConfigFromFile( std::string );

  };

}

#endif