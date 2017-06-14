#include "Linear3DChargeTaggerTypes.h"

#ifndef __CINT__
#ifndef __CLING__

namespace larlitecv {
  
  // This is how we add points to PointInfoList. We move object to avoid copy. We also tally some information.
  void PointInfoList::emplace( PointInfo&& pt ) {

    // tally
    if ( pt.planeswithcharge==(int)pt.cols.size() )
      num_pts_w_allcharge++;
    else if ( pt.planeswithcharge==0 && pt.goodpoint )
      num_pts_w_allbadch++;
    else if ( pt.planeswithcharge==0 && !pt.goodpoint )
      num_pts_w_allempty++;

    if ( pt.planeswithcharge>=2 )
      num_pts_w_majcharge++;

    if ( pt.goodpoint )
      num_pts_good++;

    // then store
    emplace_back( std::move(pt) );
  }

}

#endif
#endif

