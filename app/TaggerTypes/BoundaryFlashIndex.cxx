#include "BoundaryFlashIndex.h"

namespace larlitecv {

  BoundaryFlashIndex::BoundaryFlashIndex()
    : ivec(-1),idx(-1),popflash(NULL)
  {}
  
  BoundaryFlashIndex::BoundaryFlashIndex(int iivec, int iidx, const larlite::opflash* popf )
    : ivec(iivec),idx(iidx),popflash(popf)
  {}
  
}
