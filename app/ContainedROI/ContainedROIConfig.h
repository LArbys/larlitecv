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
		std::vector<float> min_cluster_plane_charge;
		float charge_diff_sigma;
		float charge_diff_weight;
		float time_boundary_diff_sigma;
		float time_boundary_diff_weight;
		float triarea_sigma;
		float triarea_weight;
		int max_number_rois;
		float max_roi_likelihood;
		float roi_likelihood_cutoff;
		bool generate_calib_info;
		bool draw_truth_roi;

		void setDefaults();

	};
 
	ContainedROIConfig CreateContainedROIConfig( const larcv::PSet& );

}


#endif

