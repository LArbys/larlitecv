#include "dwall.h"

#include <cmath>

namespace larlitecv {

  double dwall( const std::vector<double>& pos, int& boundary_type ) {
    std::vector<float> fpos(3);
    for (int i=0; i<3; i++)
      fpos[i] = (float)pos[i];
    return dwall( fpos, boundary_type);
  }
  
  float dwall( const std::vector<float>& pos, int& boundary_type ) {

    float dx1 = fabs(pos[0]);
    float dx2 = fabs(258-pos[0]);
    float dy1 = fabs(117.0-pos[1]);
    float dy2 = fabs(-117.0-pos[1]);
    float dz1 = fabs(pos[2]);
    float dz2 = fabs(1036.0-pos[2]);

    float dwall = 1.0e9;

    if ( dy1<dwall ) {
      dwall = dy1;
      boundary_type = 0; // top
    }
    if ( dy2<dwall ) {
      dwall = dy2;
      boundary_type = 1; // bottom
    }
    if ( dz1<dwall ) {
      dwall = dz1;
      boundary_type = 2; // upstream
    }
    if ( dz2<dwall ) {
      dwall = dz2;
      boundary_type = 3; // downstream
    }
    if ( dx1<dwall ) {
      dwall = dx1;
      boundary_type = 4; // anode
    }
    if ( dx2<dwall ) {
      dwall = dx2;
      boundary_type = 5; // cathode
    }

    return dwall;
    
  }
  
}
