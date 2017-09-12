#ifndef __crossingPointsAnaMethods_h__
#define __crossingPointsAnaMethods_h__

#include <vector>

// larlite
#include "DataFormat/opflash.h"

// larlitecv
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"

class TTree;

namespace larcv {
  class EventPixel2D;
  class Image2D;
  class ImageMeta;
}

namespace larlite {
  class trigger;
  class event_mctrack;
  class mcstep;
  class mctrack;
}


namespace larlitecv {

  class SpaceChargeMicroBooNE;

  struct TruthCrossingPointAna_t {
    // data for truth crossing points
    int start_or_end;
    int type;
    std::vector< int > imgcoord; // (row, plane cols)
    std::vector< float > crossingpt_det;    // (x,y,z) detector coordinates
    std::vector< float > crossingpt_detsce; // (x,y,z) detector coordinates w/ space charge correction
    std::vector< float > crossingpt_detsce_tyz;  // (tick, y, z) detector coordinates w/ space charge correction
    std::vector< float > crossingpt_detsce_tyz2; // (tick, y, z) detector coordinates w/ space charge correction. further inside det.
    int nplanes_w_charge;
    int flashindex;
    int mctrack_index;    
    int matched;
    int matched_type;
    int matched_recoindex;
    float matched_dist;
    void clear();
  };

  struct RecoCrossingPointAna_t {
    // data about truth matching for reco boundary points
    int type;
    int reco_flashindex; // points to position in truthcrossingptinfo_v
    int truthmatch;    
    int truthmatch_index; // nearest index
    int truthmatch_type;
    int truthmatch_flashindex;
    float truthmatch_dist;
    std::vector<float> truthmatch_detsce_tyz;
    void clear();
  };
  
  
  struct CrossingPointAnaData_t {
    int tot_proposed_crossingpoints; //< total crossing points
    int proposed_crossingpoints[7];  //< total proposed for each type [top,bot,up,down,anode,cathod,imgends]

    int tot_true_crossingpoints;
    int true_crossingpoints[7];

    int tot_flashmatched_true_crossingpoints;
    int flashmatched_true_crossingpoints[7];

    int tot_matched_crossingpoints;
    int matched_crossingpoints[7];

    int true_intime_stopmu;
    int true_intime_thrumu;

    std::vector< TruthCrossingPointAna_t > truthcrossingptinfo_v; // info stored for each true crossing point
    std::vector< RecoCrossingPointAna_t >  recocrossingptinfo_v;

    // MC track end points
    std::vector< std::vector<int> > mctrack_imgendpoint_indices; // index is mc track index. -1 or size 0 means there are no image end points
    
    bool saved_mc;
    
    CrossingPointAnaData_t() {
      clear();
    };

    void clear() {
      tot_proposed_crossingpoints = 0;
      tot_true_crossingpoints = 0;
      tot_matched_crossingpoints = 0;
      tot_flashmatched_true_crossingpoints = 0;
      for (int i=0; i<larlitecv::kNumEndTypes; i++){
	proposed_crossingpoints[i] = 0;
	true_crossingpoints[i] = 0;
	matched_crossingpoints[i] = 0;
	flashmatched_true_crossingpoints[i] = 0;
      };
      proposed_crossingpoints[6] = 0;
      true_intime_stopmu = 0;
      true_intime_thrumu = 0;
      truthcrossingptinfo_v.clear();
      recocrossingptinfo_v.clear();
      saved_mc = false;
    };

    void bindToTree( TTree* tree );
  };//end of struct

  void analyzeCrossingPoints( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta, const std::vector<larcv::Image2D>& img_v,
			      const larcv::EventPixel2D* ev_spacepoints[], const std::vector<larlite::event_opflash*>& opflash_v,
			      const larlite::trigger* ev_trigger=nullptr, const larlite::event_mctrack* ev_mctrack=nullptr );
			      
  void analyzeCrossingMCTracks( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta, const std::vector<larcv::Image2D>& img_v,
				const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack, const std::vector<larlite::event_opflash*>& opflash_v, bool printAnodeEnds );

  void analyzeCrossingDataOnly( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*>& ev_spacepoints );  

  void analyzeCrossingMatches( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*> ev_spacepoints,
			       const larcv::ImageMeta& meta, const float fMatchRadius );

  void analyzeCrossingMatches( CrossingPointAnaData_t& data,
			       const std::vector< const std::vector<larlitecv::BoundarySpacePoint>* > ev_spacepoints,
			       const larcv::ImageMeta& meta, const float fMatchRadius );

  float getTick( const std::vector<float>& step, const float trig_time=4050.0, const larlitecv::SpaceChargeMicroBooNE* psce=NULL );
  
  float getTick( const larlite::mcstep& step, const float trig_time=4050.0, const larlitecv::SpaceChargeMicroBooNE* psce=NULL );

  int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce );

  std::vector<int> getImageBoundaryCrossingPoint( const larlite::mctrack& track, std::vector<float>& crossingpt, const larcv::ImageMeta& meta,
						  const float boundary_tick_buffer, const float trig_time, const larlitecv::SpaceChargeMicroBooNE* psce );

  std::vector<float> getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time,
						 const bool startAtstart, const float max_step_size, const float fv_border, const larlitecv::SpaceChargeMicroBooNE* psce );
  
}

#endif
