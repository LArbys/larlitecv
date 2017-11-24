#ifndef __TAGGER_FLASH_MATCH_TYPES_H__
#define __TAGGER_FLASH_MATCH_TYPES_H__

#include <vector>

// LArCV
#include "DataFormat/Pixel2DCluster.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/Image2D.h"

// LArLite
#include "DataFormat/track.h"
#include "DataFormat/opflash.h"

namespace larlitecv {

  class TaggerFlashMatchData {
  public:
    typedef enum { kThruMu=0, kStopMu, kUntagged } ClusterType_t;

  TaggerFlashMatchData( ClusterType_t type, const std::vector<larcv::Pixel2DCluster>& pixels, const larlite::track& track ) 
    : m_type(type), m_pixels(pixels), m_track3d(track), m_pstart_flash(NULL), m_pend_flash(NULL), has_startflash(false), has_endflash(false), mctrackid(-1) {};
    virtual ~TaggerFlashMatchData() {};

    ClusterType_t m_type;
    std::vector<larcv::Pixel2DCluster> m_pixels;
    larlite::track m_track3d;
    const larlite::opflash* m_pstart_flash;
    const larlite::opflash* m_pend_flash;
    bool has_startflash;
    bool has_endflash;  

    larcv::ROI MakeROI( const std::vector<larcv::Image2D>& img_v, const float bbox_pad_cm=0.0, const bool iscroi_candidate=false ) const;
    void setStartFlash( const larlite::opflash* pflash ) {
      m_pstart_flash=pflash;
      if ( m_pstart_flash!=NULL )
	has_startflash=true;
      else
	has_startflash=false;
    };
    void setEndFlash( const larlite::opflash* pflash ) {
      m_pend_flash=pflash;
      if ( m_pend_flash!=NULL )
	has_endflash=true;
      else
	has_endflash=false;
    };
    bool hasStartFlash() const { return has_startflash; };
    bool hasEndFlash() const { return has_endflash; };

    int mctrackid; //< used for truth-based studies
    
  };
 
}

#endif
