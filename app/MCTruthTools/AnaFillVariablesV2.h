#ifndef __AnaFillVariablesV2_h__
#define __AnaFillVariablesV2_h__

/**********************************************************
 * AnaFillVariablesV2
 * 
 * Methods to fill Version 2 selection variables
 *
 * revision history
 * initial draft: 11/20/2017 taritree (twongj01@tufts.edu)
 *
 **********************************************************/

// stdlib
#include <vector>

// ROOT
class TTree;

// LArCV
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"

// larlite
#include "DataFormat/user_info.h"

namespace larlitecv {

  class AnaFillVariablesV2 {
  public:

    AnaFillVariablesV2() {};
    virtual ~AnaFillVariablesV2() {};
    

    virtual void bindEventTree( TTree* tree );
    virtual void bindTrackTree( TTree* tree );

    
    TTree* m_ev_tree;          //< event-level tree. assuming we do not own this object
    int   ev_croipixelarea[4]; //< event-level: pixels covered by the collection of ROIs 
    float ev_maxpe;            //< event-level: maxpe of in-time flash
    float ev_fwhm;             //< event-level: FWHM of in-time flash
    float ev_meanz;            //< event-level: Mean z position of in-time flash
    float ev_highestnufrac;    //< event-level: highest neutrino pixel fraction of all tracks

    TTree* m_track_tree;       //< track-level tree. assuming we do not own this object
    float track_maxpe_hypo;    //< track-level: maxpe of hypothesis flash
    float track_zfracdiff;     //< track-level: zfracdiff of track compared to intime flash
    float track_chi2;          //< track-level: best chi2 value to out-of-time flashes
    float track_dwall;         //< track-level: smallest dwall used to measure containment
    float track_nufrac;        //< track-level: fraction of nu pixels we've covered
    int   track_highestnufrac; //< track-level: 1 if track has highest neutrino pixel fraction in event, 0 otherwise

    void clearEventVars();
    void clearTrackVars();

    // Fill Methods
    void fillEventInfo( const larcv::EventImage2D& ev_img, const larcv::EventROI& ev_roi,
			larlite::event_user* ev_user_info, const larcv::EventImage2D* ev_segment,
			larcv::EventPixel2D* ev_allpixels_v );
    
    std::vector<larcv::Image2D> _area_counter_v;
    void fillCROIPixArea( const larcv::EventImage2D& ev_img, const larcv::EventROI& ev_roi );

    void fillIntimeFlashInfo( larlite::event_user* ev_user_info );
    void fillTrackInfo( larlite::event_user* ev_user_info,
			const larcv::EventImage2D* ev_segment=NULL,
			larcv::EventPixel2D* ev_allpixels_v=NULL );
    
    float fillTrackNuFraction( const larcv::EventImage2D* ev_segment,
			       const int itrack, larcv::EventPixel2D* ev_allpixels_v );
    
    
    
  };


}


#endif
