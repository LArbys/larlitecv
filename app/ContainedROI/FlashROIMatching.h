#ifndef __FLASHROISELECTION__
#define __FLASHROISELECTION__

/**
	// Class to Filter Out untagged ROIs using flash information

	*/

#include <vector>
#include <set>

#include "TTree.h"
#include "TFile.h"

// larcv
#include "DataFormat/ROI.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Pixel2DCluster.h"
#include "PMTWeights/PMTWireWeights.h"

// larlite
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h" // SelectionTool
#include "FhiclLite/PSet.h"

namespace larlitecv {
	
	class FlashROIMatchingConfig {
	public:
			FlashROIMatchingConfig( fcllite::PSet ps_ );
			virtual ~FlashROIMatchingConfig() {};

			void setDefaults();

			std::vector<int> beam_tick_range;
			float us_per_tick;
			float pmtflash_thresh;
			bool store_calib_data;
			std::vector<float> gain_correction;
			int min_qcluster_size;
			float pixel_threshold;
			float flash_front_boundary;
			float flash_back_boundary;
			fcllite::PSet ps;
			int verbosity;
			int qneighborhood;
			int trigger_tick;
	};

	FlashROIMatchingConfig MakeFlashROIMatchingConfigFromFile( std::string fname );

	class FlashROIMatching {

	public:
		FlashROIMatching( const FlashROIMatchingConfig& config);

		virtual ~FlashROIMatching() {};

		FlashROIMatchingConfig m_config;
		larcv::pmtweights::PMTWireWeights m_pmtweights;


		// primary routine
		std::vector< larcv::ROI > SelectFlashConsistentROIs( const std::vector<larlite::event_opflash*>& opflashes_v, const std::vector<larcv::Image2D>& img_v, 
			const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters, const std::vector< larcv::ROI >& untagged_rois,
			larcv::EventPixel2D& thrumu_clusters,
			larcv::EventPixel2D& stopmu_clusters );

		// supporting routines
    std::vector<larlite::opflash> SelectInTimeFlashes( const std::vector<larlite::event_opflash*>& opflashes_v );
    void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );
    flashana::QCluster_t MakeTPCFlashObject( const larcv::Pixel2DCluster& pixels, const larcv::Image2D& img, larcv::Pixel2DCluster& expanded_cluster );
    bool ClusterWithinRange( const std::vector<float>& wire_range, const larcv::Pixel2DCluster& pixels, const larcv::Image2D& img );

    // verbosity
    int m_verbosity;
    void SetVerbosity( int verbosity ) { m_verbosity = verbosity; }

    // flash matching algorithm manager
    flashana::FlashMatchManager m_flash_matcher;
    std::vector<flashana::FlashMatch_t> m_results;
    std::vector<flashana::FlashMatch_t> m_results_unfitted;

    // flash match data for calibration
    TFile* m_file;
    TTree* m_tree;
    int m_run;
    int m_subrun;
    int m_event;
    int m_entry;
    int m_nuflag;
    int m_tagflag;
    int m_passcuts;
    int m_clustertype;
    float m_totalpe;
    float m_rawtotpe;
    float m_predictpe;
    float m_score;
    float m_fraction_nupixels;
    float m_flashchi2;
    float m_flash_range[2];
    float m_flash_hypothesis[32];
    float m_raw_hypothesis[32];
    float m_measured[32];
    void writeCalibTree() { if ( m_config.store_calib_data) { m_file->cd(); m_tree->Write(); } }

    larcv::Image2D* ptr_segimg;
    void setNeutrinoYPlaneSegmentImage( larcv::Image2D* seg ) { ptr_segimg = seg; }
    int NumOfSegPixels( const larcv::Image2D& img );
		int NumOfSegPixelOverlap( const larcv::Pixel2DCluster& pixels, const larcv::ImageMeta& meta );
		void setRSE( int entry, int run, int subrun, int event ) { m_entry = entry; m_run = run; m_subrun = subrun; m_event = event; }

		struct CandidateFlashMatchedROI_t {
			// struct ties together TPC flash objects to the pixel2d clusters from which they derive.
			// they also track their compatibilty to the in-time beam flash(es)
			std::vector<int> flash_idxs; // indices of flash objects the 'track object' was compatible with
			int qcluster_idx; // index of QCluster_t object in flashmatchmanager
			int source_type; // 0: untagged clusters,  1: thrumu clusters, 2: stopmu clusters
			int source_idx;  // position in the source's container
			flashana::QCluster_t qcluster;
			std::vector< larcv::Pixel2DCluster > expanded_clusters; // pixels in each plane. found by expanding search for charge around base track objects
		};

		std::vector<CandidateFlashMatchedROI_t> GenerateFlashTPCObjects( const std::vector<larlite::opflash>& beam_flashes, 
			const std::vector< std::vector<float> >& flash_wire_range,
			const std::vector< std::vector<larcv::Pixel2DCluster> >& untagged_clusters, const std::vector< larcv::ROI >& untagged_rois,
			larcv::EventPixel2D& thrumu_clusters,
			larcv::EventPixel2D& stopmu_clusters,
			const std::vector<larcv::Image2D>&  img_v );


		std::vector<flashana::FlashMatch_t> GenerateFittedFlashHypotheses( std::vector< CandidateFlashMatchedROI_t> & candidates );
		std::vector<flashana::FlashMatch_t> GenerateUnfittedFlashHypotheses( std::vector< CandidateFlashMatchedROI_t> & candidates );		

		void SelectCandidateROIs( std::vector<CandidateFlashMatchedROI_t>& candidates, 
		std::vector< flashana::FlashMatch_t>& matches, std::set< std::pair<int,int> >& cluster_indices );

		larcv::ROI GenerateClusterROI( const std::vector<larcv::Pixel2DCluster>& cluster, const std::vector<larcv::Image2D>& img_v );

	
	};

}

#endif
