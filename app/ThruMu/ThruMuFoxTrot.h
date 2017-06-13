#ifndef __THRUMU_FOXTROT_H__
#define __THRUMU_FOXTROT_H__

/* ============================================================
 * 
 * ThruMuFoxTrot
 *
 * This tracker uses the fox trot stepping algorithm.
 * We attempt to connect two points under the through-going muon assumption.  
 * We make a custom FoxTrot instance that chooses points based on a straight track 
 * after space charge corrections.
 * 
 * We also provide flash information for anode/cathode ends, this way we can evaluate for poor
 * tracks. (maybe move this to another system)
 *
 * ===========================================================*/

#include <vector>

// larlite
#include "BasicTool/GeoAlgo/GeoAlgo.h"

// larcv
#include "DataFormat/Image2D.h" 

#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "SCE/ReverseSCE.h"
#include "SCE/SpaceChargeMicroBooNE.h"
#include "ChargeSegmentAlgos/FoxTrotLead.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgo.h"
#include "T3DMerge/T3DCluster.h"
#include "T3DMerge/T3DPCMerge.h"

#include "ThruMuFoxTrotConfig.h"


namespace larlitecv {

  class ThruMuFoxTrotLead : public FoxTrotLead {
  public:
    ThruMuFoxTrotLead() {};
    virtual ~ThruMuFoxTrotLead() {};

    void setEndPoints( const std::vector<double>& start, const std::vector<double>& end, const double shiftx=0 );
    void setEndPoints( const std::vector<float>& start, const std::vector<float>& end, const float shiftx=0 );    

    std::vector<double> getSCECorrectedPos( const std::vector<double>& pos ) { return m_sce.getOriginalPos( pos ); };
    std::vector<double> getSCEAppliedPos( const std::vector<double>& pos );
    double getDist2RealLine( const std::vector<double>& pt );
    std::vector<double> getClosestPointOnRealLine( const std::vector<double>& pt );
    const std::vector<double>& getRealLineDir() const;
      
  protected:

    bool _chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past, const FoxTrotTrackerAlgoConfig& config,
			      const std::vector<larcv::Image2D>& stepped_v, int& best_seg_idx );
    
    //float dist2realline( std::vector<double>& image_pt );
    
    // we'll use the larlite geoalgo functions
    geoalgo::Line image_line3d; ///< line according to the image
    geoalgo::Line real_line3d;  ///< line in real-space, after spacecharge correction
    std::vector<double> m_realdir;
    double m_shiftx;


    larlitecv::ReverseSCE m_sce; // get the offsets
    larlitecv::SpaceChargeMicroBooNE m_sce_forward;
    geoalgo::GeoAlgo m_geoalgo;
  };

  class ThruMuFoxTrot {
  public:
    ThruMuFoxTrot( const ThruMuFoxTrotConfig& cfg );
    virtual ~ThruMuFoxTrot() {};

    BMTrackCluster3D findThruMuTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );
    
    FoxTrack findThruMuFoxTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );

    const std::vector<int>& getStepTypes() { return m_step_types; };
    
  protected:
    
    ThruMuFoxTrotConfig m_config;
    FoxTrotTrackerAlgo m_tracker;
    ThruMuFoxTrotLead m_lead;
    std::vector<int> m_step_types; // 0=trot tracker 1=algo

    FoxStep makeHypothesisStep( const FoxStep& current_step, const double step_size );    
    bool doesHypothesisStepSeeCharge( FoxStep& hypostep, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );
    double distanceToTheEnd( FoxStep& current, const BoundarySpacePoint& endpoint );

    // yet to be implemented
    void cleanTrackEnd( FoxTrack& track, float cos_requirement );
    std::vector< T3DCluster > makeRealAndImagePaths( FoxTrack& trot );
    std::vector< std::vector<T3DCluster> > breakPathsIntoStraightSegments( std::vector<T3DCluster>& paths );
    std::vector< std::vector<T3DCluster> > m_track_splits;
    T3DPCMerge m_pcmerge;
    
    
    std::vector< T3DCluster > m_real_paths;  //< paths of segments that have been sce-corrected, and broken into straight line segments
    std::vector< T3DCluster > m_image_paths; //< paths of segments based on points in the images, points correspond to partner segments in m_real_paths
    geoalgo::GeoAlgo m_geoalgo;    
    
  };
  
  
}

#endif
