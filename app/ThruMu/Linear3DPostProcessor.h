#ifndef __LINEAR3DPOSTPROCESSOR_H__
#define __LINEAR3DPOSTPROCESSOR_H__

/*! \brief Post-process BMTrackCluster3D objects from Linear3DFitter
 *
 *
 * This class processes a vector of BMTrackCluster3D objects.
 * It tries to find tracks that are a subsegment of a larger track.
 * Also, if it can merge two tracks, it will try to do that too.
 *
 * Basically, it tries to clean up the output and make up for failures of 
 * the Linear3DFitter.
 *
 * author(s): Taritree Wongjirad (taritree@mit.edu)
 *
 * revisions
 * 2/14/2017: first writing
 */

#include <vector>

// larlitecv/app/Thrumu
#include "BMTrackCluster3D.h"

namespace larlitecv {

  class Linear3DPostProcessor {

  public:
    Linear3DPostProcessor() {};
    virtual ~Linear3DPostProcessor() {};

    std::vector< BMTrackCluster3D > process( std::vector< BMTrackCluster3D >& tracks_v );
    

  };

}

#endif
