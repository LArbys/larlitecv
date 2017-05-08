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
#include "ChargeSegmentAlgos/FoxTrotLead.h"
#include "ChargeSegmentAlgos/FoxTrotTrackerAlgo.h"

#include "ThruMuFoxTrotConfig.h"


namespace larlitecv {

  class ThruMuFoxTrotLead : public FoxTrotLead {
  public:
    ThruMuFoxTrotLead() {};
    virtual ~ThruMuFoxTrotLead() {};

    void setEndPoints( const std::vector<double>& start, const std::vector<double>& end, const double shiftx=0 );
    void setEndPoints( const std::vector<float>& start, const std::vector<float>& end, const float shiftx=0 );    
    
  protected:

    bool _chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past, const FoxTrotTrackerAlgoConfig& config, int& best_seg_idx );
    
    //float dist2realline( std::vector<double>& image_pt );
    
    // we'll use the larlite geoalgo functions
    geoalgo::Line image_line3d; ///< line according to the image
    geoalgo::Line real_line3d;  ///< line in real-space, after spacecharge correction
    double m_shiftx;


    ReverseSCE m_sce; // get the offsets
    geoalgo::GeoAlgo m_geoalgo;
  };

  class ThruMuFoxTrot {
  public:
    ThruMuFoxTrot( const ThruMuFoxTrotConfig& cfg );
    virtual ~ThruMuFoxTrot() {};

    BMTrackCluster3D findThruMuTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );
    
    FoxTrack findThruMuFoxTrack( const BoundarySpacePoint& pt_a, const BoundarySpacePoint& pt_b, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );
			      
    
  protected:
    
    ThruMuFoxTrotConfig m_config;
    FoxTrotTrackerAlgo m_tracker;
    ThruMuFoxTrotLead m_lead;
    
  };
  
  
}

#endif
