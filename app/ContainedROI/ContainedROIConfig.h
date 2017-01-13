#ifndef __CONTAINEDROI_CONFIG__
#define __CONTAINEDROI_CONFIG__

// larcv
#include "Base/PSet.h"

namespace larlitecv {

	class ContainedROIConfig {
	public:
		ContainedROIConfig() {};
		virtual ~ContainedROIConfig() {};

		float pixel_threshold;
		int min_cluster_size;

		void setDefaults();

	};
 
	ContainedROIConfig CreateContainedROIConfig( const larcv::PSet& );

}


#endif

