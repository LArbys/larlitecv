#ifndef __TAGGER_CROI_ALGO_CONFIG_H__
#define __TAGGER_CROI_ALGO_CONFIG_H__

// larcv
#include "Base/PSet.h"

namespace larlitecv {

  class TaggerCROIAlgoConfig {
  public:
    TaggerCROIAlgoConfig() {};
    virtual ~TaggerCROIAlgoConfig() {};

    larcv::PSet sidetagger_pset;
    larcv::PSet flashtagger_pset;
    
  };

}

#endif