#ifndef __TAGGER_CROI_ALGO_CONFIG_H__
#define __TAGGER_CROI_ALGO_CONFIG_H__

// larcv
#include "Base/PSet.h"

// larlitecv
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "StopMu/StopMuFilterSpacePoints.h"
#include "StopMu/StopMuClusterConfig.h"

namespace larlitecv {

  class TaggerCROIAlgoConfig {
  public:
    TaggerCROIAlgoConfig() {};
    virtual ~TaggerCROIAlgoConfig() {};

    larlitecv::ConfigBoundaryMuonTaggerAlgo  sidetagger_cfg;
		larlitecv::FlashMuonTaggerConfig         flashtagger_cfg;
		larlitecv::StopMuFilterSpacePointsConfig stopmu_filterpts_cfg;
		larlitecv::StopMuClusterConfig           stopmu_cluster_cfg;

    static TaggerCROIAlgoConfig makeConfigFromFile( std::string );

  };

}

#endif