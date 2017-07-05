#include "TaggerFlashMatchAlgoConfig.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larlitecv {

	const std::string TaggerFlashMatchAlgoConfig::m_flashman_default = "FlashMatchManager: {AllowReuseFlash: true CustomAlgo: [0,0] FlashFilterAlgo: \"\" HypothesisAlgo: \"PhotonLibHypothesis\" MatchAlgo: \"QLLMatch\" ProhibitAlgo: \"\" StoreFullResult: true TPCFilterAlgo: \"TimeCompatMatch\" Verbosity: 02 }";

	TaggerFlashMatchAlgoConfig::TaggerFlashMatchAlgoConfig():
	 m_flashmatch_config(fcllite::PSet("FlashMatchManager",m_flashman_default)) {
		verbosity = 0;
		qcluster_stepsize = 0.3;
		MeV_per_cm = 2.3;
		fudge_factor = 33333.0;
		gain_correction.resize(32,1.0);
		pmtflash_thresh = 3;
		use_gaus2d = false;
	}

	TaggerFlashMatchAlgoConfig TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( const larcv::PSet& pset ) {
		TaggerFlashMatchAlgoConfig cfg;
		cfg.verbosity         = pset.get<int>("Verbosity");
		cfg.qcluster_stepsize = pset.get<float>("QClusterStepSize");
		cfg.MeV_per_cm        = pset.get<float>("MeV_per_cm");
		cfg.fudge_factor      = pset.get<float>("FudgeFactor");
		cfg.pmtflash_thresh   = pset.get<float>("PMTFlashThreshold");
		cfg.use_gaus2d        = pset.get<bool>("UseGaus2D");
		fcllite::PSet vox( "Manager", pset.data_string() );
		cfg.m_flashmatch_config = vox.get<fcllite::PSet>("TaggerFlashMatchAlgo");
		return cfg;
	}
}
