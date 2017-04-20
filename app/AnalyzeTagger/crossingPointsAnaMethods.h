#ifndef __crossingPointsAnaMethods_h__
#define __crossingPointsAnaMethods_h__

#include <vector>
#include "DataFormat/opflash.h"

class TTree;

namespace larcv {
  class EventPixel2D;
  class ImageMeta;
}

namespace larlite {
  class trigger;
  class event_mctrack;
}

namespace larlitecv {

  struct CrossingPointAnaData_t {
    int tot_proposed_crossingpoints; //< total crossing points
    int proposed_crossingpoints[7];  //< total proposed for each type [top,bot,up,down,anode,cathod,imgends]

    int tot_true_crossingpoints;
    int true_crossingpoints[6];

    int tot_flashmatched_true_crossingpoints;
    int flashmatched_true_crossingpoints[6];

    int tot_matched_crossingpoints;
    int matched_crossingpoints[6];

    int true_intime_stopmu;
    int true_intime_thrumu;

    std::vector< int > start_type;
    std::vector< int > end_type;
    std::vector< std::vector<int> > start_pixels;
    std::vector< std::vector<float> > start_crossingpts;
    std::vector< std::vector<int> > end_pixels;
    std::vector< std::vector<float> > end_crossingpts;
    bool saved_mc;
    
    CrossingPointAnaData_t() {
      clear();
    };

    void clear() {
      tot_proposed_crossingpoints = 0;
      tot_true_crossingpoints = 0;
      tot_matched_crossingpoints = 0;
      tot_flashmatched_true_crossingpoints = 0;
      for (int i=0; i<6; i++){
	proposed_crossingpoints[i] = 0;
	true_crossingpoints[i] = 0;
	matched_crossingpoints[i] = 0;
	flashmatched_true_crossingpoints[i] = 0;
      };
      proposed_crossingpoints[6] = 0;
      true_intime_stopmu = 0;
      true_intime_thrumu = 0;
      start_type.clear();
      end_type.clear();
      start_pixels.clear();
      start_crossingpts.clear();
      end_pixels.clear();
      end_crossingpts.clear();
      saved_mc = false;
    };

    void bindToTree( TTree* tree );
  };//end of struct

  void analyzeCrossingPoints( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta, const larcv::EventPixel2D* ev_spacepoints[], const std::vector<larlite::event_opflash*>& opflash_v,
			      const larlite::trigger* ev_trigger=nullptr, const larlite::event_mctrack* ev_mctrack=nullptr );
			      
  void analyzeCrossingMCTracks( CrossingPointAnaData_t& data, const larcv::ImageMeta& meta,
				const larlite::trigger* ev_trigger, const larlite::event_mctrack* ev_mctrack, const std::vector<larlite::event_opflash*>& opflash_v, bool printAnodeEnds );

  void analyzeCrossingDataOnly( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*>& ev_spacepoints );  

  void analyzeCrossingMatches( CrossingPointAnaData_t& data, std::vector<larcv::EventPixel2D*> ev_spacepoints, const larcv::ImageMeta& meta );  
  
}

#endif
