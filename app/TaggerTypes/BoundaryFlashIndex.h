#ifndef __BOUNDARY_FLASH_INDEX_H__
#define __BOUNDARY_FLASH_INDEX_H__

/* =====================================================
 * BoundaryFlashIndex
 *
 * Struct for book keeping purposes that associates
 *  BoundarySpacePoint objects to
 *  individual instances of opflash in 
 *  vector< event_opflash > 
 * =====================================================*/

#include "DataFormat/opflash.h"

namespace larlitecv {

  class BoundaryFlashIndex {
  public:
    BoundaryFlashIndex();
    BoundaryFlashIndex(int iivec, int iidx, const larlite::opflash* popf=NULL);
    virtual ~BoundaryFlashIndex() {};

    int ivec; //< index of container: typically 0=beam flashes, 1=cosmic flashes
    int idx;  //< index of opflash instance in the container
    const larlite::opflash* popflash; //< pointer to it

  };

}

#endif
