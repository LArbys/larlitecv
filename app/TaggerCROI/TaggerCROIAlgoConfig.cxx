#include "TaggerCROIAlgoConfig.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

namespace larlitecv {
	
	TaggerCROIAlgoConfig::TaggerCROIAlgoConfig() 
	 : input_write_cfg(larcv::PSet("input")), 
	 	 thrumu_write_cfg(larcv::PSet("thrumu")), 
	 	 stopmu_write_cfg(larcv::PSet("stopmu")),
	 	 croi_write_cfg(larcv::PSet("croi"))
	{

	}

	TaggerCROIAlgoConfig TaggerCROIAlgoConfig::makeConfigFromFile( std::string filepath ) {

		TaggerCROIAlgoConfig cfg;

		larcv::PSet root_pset           = larcv::CreatePSetFromFile( filepath );
		larcv::PSet tagger_pset         = root_pset.get<larcv::PSet>("TaggerCROI");
		larcv::PSet sidetagger_pset     = tagger_pset.get<larcv::PSet>("BMTSideTagger");
		larcv::PSet flashtagger_pset    = tagger_pset.get<larcv::PSet>("BMTFlashTagger");
		larcv::PSet stopmu_filter_pset  = tagger_pset.get<larcv::PSet>("StopMuSpacePointsFilter");
		larcv::PSet stopmu_cluster_pset = tagger_pset.get<larcv::PSet>("StopMuCluster");
		larcv::PSet untagged_cluster_ps = tagger_pset.get<larcv::PSet>("ContainedGroupAlgo");
		larcv::PSet untagged_match_pset = tagger_pset.get<larcv::PSet>("TaggerFlashMatchAlgo");

		cfg.input_write_cfg = tagger_pset.get<larcv::PSet>("InputWriteConfig");
		cfg.thrumu_write_cfg = tagger_pset.get<larcv::PSet>("ThruMuWriteConfig");

		cfg.sidetagger_cfg       = larlitecv::MakeConfigBoundaryMuonTaggerAlgoFromPSet( sidetagger_pset );
    cfg.flashtagger_cfg      = larlitecv::MakeFlashMuonTaggerConfigFromPSet( flashtagger_pset );
    cfg.stopmu_filterpts_cfg = larlitecv::MakeStopMuFilterSpacePointsConfigFromPSet( stopmu_filter_pset );
    cfg.stopmu_cluster_cfg   = larlitecv::makeStopMuClusterConfigFromPSet( stopmu_cluster_pset );
    cfg.untagged_cluster_cfg = larlitecv::ClusterGroupAlgoConfig::MakeClusterGroupAlgoConfigFromPSet( untagged_cluster_ps );
    cfg.croi_selection_cfg   = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( untagged_match_pset );

		return cfg;
	}

}