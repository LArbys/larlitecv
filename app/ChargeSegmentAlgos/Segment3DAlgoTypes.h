#ifndef __SEGMENT3D_ALGO_TYPES__
#define __SEGMENT3D_ALGO_TYPES__

namespace larlitecv {

  struct Segment3D_t {
    float start[3];
    float end[3];
  };

  struct Segment2D_t {
    int row_high;
    int col_high;
    int row_low;
    int col_low;
  };

}

#endif
