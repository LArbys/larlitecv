//Skeletonization
//Class written by Teo Popa
//2006
//implement based on Nicu Cornea CThinning code given by 
//Palagyi and Kuba in "A Parallel 3D 12-Subiteration CThinning Algorithm", 1999
//

//
// This code modded by vic
//


#ifndef __VTHINNING_H__
#define __VTHINNING_H__

#include <memory.h>
#include <vector>

#include "VPoint.h"

using namespace std;

#define OBJECT 1
#define D_BORDER 10
#define SIMPLE 20

namespace llcv {
  enum Direction {
    UP_SOUTH = 0, 
    NORT_EAST,
    DOWN_WEST, 
  
    SOUTH_EAST,
    UP_WEST, 
    DOWN_NORTH, 
  
    SOUTH_WEST,
    UP_NORTH, 
    DOWN_EAST, 
  
    NORT_WEST,
    UP_EAST, 
    DOWN_SOUTH,
  
    UP,
    DOWN, 
    EAST, 
    WEST, 
    NORTH, 
    SOUTH
  };

  class VThinning{

  public:
    VThinning();
    ~VThinning() {}

    int Execute(unsigned char *volin,int SizeX,int  SizeY,int  SizeZ, vector<VPoint> *pExtractedIndexes);
    bool MatchesATemplate(unsigned char n[3][3][3]);
    bool TransformNeighborhood(unsigned char n[3][3][3], char direction,unsigned char USn[3][3][3]); 
    bool markBoundaryInDirection(unsigned char *vol, int L, int M, int N,char direction);
    inline void CopyNeighborhoodInBuffer(unsigned char *vol, int L, int M, int N, 
					 int idx, unsigned char nb[3][3][3], 
					 bool changeValues = true)
    {
      int nidx;
      char i, j, k, ii;
          
      ii = 0;
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    nidx = idx + volNeighbors[ii];
        	
	    if(!changeValues) {
	      nb[i][j][k] = vol[nidx];
	    }
	    else {
	      if(vol[nidx] != 0) {
		nb[i][j][k] = OBJECT;
	      }
	      else {
		nb[i][j][k] = 0;
	      }
	    }
        	
	    ii++;
	  }
	}
      }
          
      return;
    }

  public:
    int volNeighbors[27]; 

  };
}
#endif
