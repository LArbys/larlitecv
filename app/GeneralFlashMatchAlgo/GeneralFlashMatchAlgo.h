/* Function to compare a flash of light with a track using the flash matching infrastructure. 

   This is meant to introduce flash matching into the tagger's operations. */
#ifndef __GENERAL_FLASHMATCH_ALGO_H__
#define __GENERAL_FLASHMATCH_ALGO_H__

#include <vector>
#include <cmath>

// larlite 
#include "OpT0Finder/Algorithms/QLLMatch.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/opflash.h"
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "OpT0Finder/Algorithms/LightPath.h"
#include "DataFormat/mctrack.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/TimeService.h" // For the G4 to electronics time.

// larcv
#include "Base/PSet.h"
#include "PMTWeights/PMTWireWeights.h"

// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "GeneralFlashMatchAlgoConfig.h"

namespace larlitecv {

  // define the class
  class GeneralFlashMatchAlgo {
    
    GeneralFlashMatchAlgo();
    
  public:
    GeneralFlashMatchAlgo( GeneralFlashMatchAlgoConfig& config );
    virtual ~GeneralFlashMatchAlgo() {};

    // This might need more information based on the needs of converting a flash to an opflash.
    std::vector < int > non_anodecathode_flash_idx( const std::vector < larlite::opflash* >& opflash_v, std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v );

    std::vector < larlite::opflash* > generate_single_denomination_flash_list( std::vector< larlite::opflash > full_opflash_v, std::vector < int > idx_v);

    std::vector < larlite::track >  generate_tracks_between_passes( const std::vector< larlitecv::BMTrackCluster3D > trackcluster3d_v);

    void ExpandQClusterStartingWithLarliteTrack( flashana::QCluster_t qcluster, const larlite::track& larlite_track, double extension, bool extend_start, bool extend_end ); 

    bool isNearActiveVolumeEdge( ::geoalgo::Vector pt, double d );
    
    void GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange );

    flashana::Flash_t GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster );

    larlite::opflash MakeOpFlashFromFlash( const flashana::Flash_t& Flash );

    larlite::opflash GetGaus2DPrediction( const larlite::opflash& flash );
    
    flashana::Flash_t MakeDataFlash( const larlite::opflash& flash );

    larlite::opflash make_flash_hypothesis(const larlite::track input_track);

    std::vector<larlite::opflash> make_flash_hypothesis_vector( const std::vector<larlite::track>& input_track_v );

    float generate_chi2_in_track_flash_comparison(const flashana::QCluster_t qcluster, const larlite::opflash data_opflash, int flash_prod_idx);
    
    const GeneralFlashMatchAlgoConfig& getConfig() { return m_config; };
    
    void setVerbosity( int v ) { m_verbosity = v; };

  private:
    const GeneralFlashMatchAlgoConfig m_config;   //< configuration class
    larcv::pmtweights::PMTWireWeights m_pmtweights;
    flashana::FlashMatchManager m_flash_matcher; //< Generates Flash Hypothesis producer
    //SpaceChargeMicroBooNE m_sce;           //< Space Charge Effect Calculator
    int m_verbosity;
    std::map<int,int> m_opch_from_opdet;

  };


} 

#endif

