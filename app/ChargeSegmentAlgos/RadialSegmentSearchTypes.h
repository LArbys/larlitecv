#ifndef __RADIAl_SEGMENT_SEARCH_TYPES_H__
#define __RADIAl_SEGMENT_SEARCH_TYPES_H__

#include <vector>
#include "DataFormat/Pixel2D.h"

namespace larlitecv {

  class RadialHit_t {
  public:
    RadialHit_t() {};
    virtual ~RadialHit_t() {};
    int start_idx;
    int end_idx;
    int max_idx;
    float maxval;
    int min_row;
    int max_row;
    std::vector<larcv::Pixel2D> pixlist;
    void reset() {
      start_idx = -1;
      end_idx = -1;
      max_idx = -1;
      maxval = -1;
      min_row = -1;
      max_row = -1;
      pixlist.clear();
      pixlist.reserve(10);
    };
    void update_max( int idx, float val ) {
      if ( val > maxval ) {
	max_idx = idx;
	maxval = val;
      }
    };
    void add_pixel( const larcv::Pixel2D& pix ) {
      if ( min_row<0 || min_row>(int)pix.Y() )
	min_row = pix.Y();
      if ( max_row<0 || max_row<(int)pix.Y() )
	max_row = pix.Y();
      pixlist.push_back( pix );
    };
    void set_end( int idx ) {
      end_idx = idx;
      end_idx -= start_idx;
      max_idx -= start_idx;
      //std::cout << "set end: start=" << start_idx << " max=" << max_idx << " end=" << end_idx << " pixlist.size=" << pixlist.size() << std::endl;
      if (max_idx>=(int)pixlist.size())
	max_idx = pixlist.size()-1;      
      start_idx = 0;
    };

  };
  

}

#endif
