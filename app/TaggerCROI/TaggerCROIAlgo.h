#ifndef __TAGGER_CROI_ALGO_H__
#define __TAGGER_CROI_ALGO_H__

#include "TaggerCROIAlgoConfig.h"

namespace larlitecv {

  class TaggerCROIAlgo {

    TaggerCROIAlgo() {};

  public:

  	TaggerCROIAlgo( const TaggerCROIAlgoConfig& config ) : m_config(config) {};
    virtual ~TaggerCROIAlgo() {};

    TaggerCROIAlgoConfig m_config;

  };
  

}

#endif
