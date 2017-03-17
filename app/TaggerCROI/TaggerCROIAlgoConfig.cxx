#include "TaggerCROIAlgoConfig.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

namespace larlitecv {
	
	TaggerCROIAlgoConfig TaggerCROIAlgoConfig::makeConfigFromFile( std::string filepath ) {

		TaggerCROIAlgoConfig cfg;

		larcv::PSet root_pset           = larcv::CreatePSetFromFile( filepath );
		larcv::PSet tagger_pset         = root_pset.get<larcv::PSet>("TaggerCROI");
		larcv::PSet sidetagger_pset     = tagger_pset.get<larcv::PSet>("BMTSideTagger");
		larcv::PSet flashtagger_pset    = tagger_pset.get<larcv::PSet>("BMTFlashTagger");
		larcv::PSet stopmu_filter_pset  = tagger_pset.get<larcv::PSet>("StopMuSpacePointsFilter");
		larcv::PSet stopmu_cluster_pset = tagger_pset.get<larcv::PSet>("StopMuCluster");
		
		cfg.sidetagger_cfg       = larlitecv::MakeConfigBoundaryMuonTaggerAlgoFromPSet( sidetagger_pset );
    cfg.flashtagger_cfg      = larlitecv::MakeFlashMuonTaggerConfigFromPSet( flashtagger_pset );
    cfg.stopmu_filterpts_cfg = larlitecv::MakeStopMuFilterSpacePointsConfigFromPSet( stopmu_filter_pset );
    cfg.stopmu_cluster_cfg   = larlitecv::makeStopMuClusterConfigFromPSet( stopmu_cluster_pset );

		return cfg;
	}

}