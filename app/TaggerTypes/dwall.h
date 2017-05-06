#ifndef __DWALL_H__
#define __DWALL_H__

#include <vector>

namespace larlitecv {
  
  float dwall( const std::vector<float>& pos, int& boundary_type );
  double dwall( const std::vector<double>& pos, int& boundary_type );  
  
}

#endif
