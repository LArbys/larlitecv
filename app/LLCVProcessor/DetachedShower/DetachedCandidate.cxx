#ifndef __DETACHEDCANDIDATE_CXX__
#define __DETACHEDCANDIDATE_CXX__

#include "DetachedCandidate.h"

namespace llcv {

  void DetachedCandidate::Insert(const llcv::DetachedCluster& dc,int plane) {
    if (plane >=3 or plane<0) throw std::runtime_error("Invalid plane");
    _dcluster_v[plane] = dc;
  }

  void DetachedCandidate::Move(llcv::DetachedCluster&& dc,int plane) {
    if (plane >=3 or plane<0) throw std::runtime_error("Invalid plane");
    _dcluster_v[plane] = std::move(dc);
  }

}
#endif
