#include "T3D2LarliteTrack.h"

namespace larlitecv {

  larlite::track T3D2LarliteTrack( const T3DCluster& track ) {
    larlite::track lltrack;
    std::cout << "convert track with " << track.getPath().size() << " and " << track.getPathDir().size() << " points into larlite" << std::endl;
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
	std::cout << "pathdir size=" << track.getPathDir().size() << std::endl;
	std::cout << "WTF" << std::endl;
	for (int i=0; i<3; i++)
	  dir[i] = track.getPathDir().back()[i];
      }
      lltrack.add_vertex( pos );
      lltrack.add_direction( dir );
    }
    return lltrack;
  }
  
}
