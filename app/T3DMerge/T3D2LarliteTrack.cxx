#include "T3D2LarliteTrack.h"

namespace larlitecv {

  larlite::track T3D2LarliteTrack( const T3DCluster& track ) {
    larlite::track lltrack;
    for ( int ipt=0; ipt<(int)track.getPath().size(); ipt++ ) {
      TVector3 pos;
      for (int i=0; i<3; i++)
	pos[i] = track.getPath()[ipt][i];
      TVector3 dir;
      if ( ipt<(int)track.getPath().size()-1 ) {
	for (int i=0; i<3; i++)
	  dir[i] = track.getPathDir()[ipt][i];
      }
      else {
	for (int i=0; i<3; i++)
	  dir[i] = track.getPathDir().back()[i];
      }
      lltrack.add_vertex( pos );
      lltrack.add_direction( dir );
    }
    return lltrack;
  }
  
}
