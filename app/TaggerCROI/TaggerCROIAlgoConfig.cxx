#include "TaggerCROIAlgoConfig.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

namespace larlitecv {
	
	TaggerCROIAlgoConfig TaggerCROIAlgoConfig::makeConfigFromFile( std::string filepath ) {

		TaggerCROIAlgoConfig cfg;

		larcv::PSet root_pset   = larcv::CreatePSetFromFile( filepath );
		larcv::PSet tagger_pset = root_pset.get<larcv::PSet>("TaggerCROI");
		cfg.sidetagger_pset     = tagger_pset.get<larcv::PSet>("BMTSideTagger");
		cfg.flashtagger_pset    = tagger_pset.get<larcv::PSet>("BMTFlashTagger");

		return cfg;
	}

}