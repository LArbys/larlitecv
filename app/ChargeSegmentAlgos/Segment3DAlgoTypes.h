#ifndef __SEGMENT3D_ALGO_TYPES__
#define __SEGMENT3D_ALGO_TYPES__

#include <vector>

namespace larlitecv {

  struct Segment3D_t {
    float start[3];
    float end[3];
    std::vector<float> plane_frac_w_charge;
    Segment3D_t() {
      plane_frac_w_charge.resize(3,0.0);
    };
  };

  struct Segment2D_t {
    int row_high;
    int col_high;
    int row_low;
    int col_low;
    float frac_w_charge;
    int npix_w_charge;
    int npix_tot;
    Segment2D_t() {
      frac_w_charge = 0.0;
      npix_w_charge = 0;
      npix_tot = 0;
    };
  };

}

#endif
