#ifndef __T3D_2_LARLITE_TRACK_H__
#define __T3D_2_LARLITE_TRACK_H__

#include "DataFormat/track.h"

#include "T3DCluster.h"

namespace larlitecv {

  larlite::track T3D2LarliteTrack( const T3DCluster& track );

}

#endif
