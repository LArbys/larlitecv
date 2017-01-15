#include "ContainedROIConfig.h"

namespace larlitecv {

	void ContainedROIConfig::setDefaults() {
		pixel_threshold = 10.;
		min_cluster_size = 5;
		min_cluster_plane_charge = std::vector<float>(3,1000.0);
		charge_diff_sigma = 0.2;
		charge_diff_weight = 1.0;
		time_boundary_diff_sigma = 10.;
		time_boundary_diff_weight = 1.0;
		triarea_sigma = 3.0;
		triarea_weight = 0.5;
		max_number_rois = 5;
		max_roi_likelihood = 10000.0;
		roi_likelihood_cutoff = 100.0;
		generate_calib_info = false;
	}

	ContainedROIConfig CreateContainedROIConfig( const larcv::PSet& ps ) {

		ContainedROIConfig config;
		config.pixel_threshold   = ps.get<float>("PixelThreshold");
		config.min_cluster_size  = ps.get<int>("MinClusterSize");
		config.min_cluster_plane_charge = ps.get< std::vector<float> >("MinClusterPlaneCharge");
		config.charge_diff_sigma = ps.get<float>("ChargeDiffSigma");
		config.charge_diff_weight = ps.get<float>("ChargeDiffWeight");
		config.time_boundary_diff_sigma = ps.get<float>("TimeBoundaryDiffSigma");
		config.time_boundary_diff_weight = ps.get<float>("TimeBoundaryDiffWeight");
		config.triarea_sigma = ps.get<float>("TriAreaSigma");
		config.triarea_weight = ps.get<float>("TriAreaWeight");
		config.max_number_rois = ps.get<int>("MaxNumberOfROIs");
		config.max_roi_likelihood = ps.get<float>("MaxROIlikelihood");
		config.roi_likelihood_cutoff = ps.get<float>("ROIlikelihoodCutoff");
		config.generate_calib_info = ps.get<bool>("GenerateCalibInfo",false);

		return config;
	}


}