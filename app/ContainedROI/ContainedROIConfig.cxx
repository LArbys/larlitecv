#include "ContainedROIConfig.h"

namespace larlitecv {

	void ContainedROIConfig::setDefaults() {
		pixel_threshold = 10.;
		min_cluster_size = 5;
	}

	ContainedROIConfig CreateContainedROIConfig( const larcv::PSet& ps ) {
		ContainedROIConfig config;
		config.pixel_threshold   = ps.get<float>("PixelThreshold");
		config.min_cluster_size  = ps.get<int>("MinClusterSize");

		return config;
	}


}