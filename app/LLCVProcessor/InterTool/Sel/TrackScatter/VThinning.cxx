#ifndef __VTHINNING_CXX__
#define __VTHINNING_CXX__

#include "VThinning.h"

#include <iostream>


namespace llcv {

  VThinning::VThinning() {
    for(int i=0;i<27;i++)
      volNeighbors[i]=0;
  }


  bool VThinning::MatchesATemplate(unsigned char n[3][3][3]) {
    // T1
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT))
       &&
       ((n[0][0][2] + n[1][0][2] + n[2][0][2] + 
	 n[0][1][2] + n[1][1][2] + n[2][1][2] + 
	 n[0][2][2] + n[1][2][2] + n[2][2][2]) == 0)
       &&
       ((n[0][0][0] + n[1][0][0] + n[2][0][0] + n[0][1][0] + n[2][1][0] + n[0][2][0] + n[1][2][0] + n[2][2][0] + 
	 n[0][0][1] + n[1][0][1] + n[2][0][1] + n[0][1][1] + n[2][1][1] + n[0][2][1] + n[1][2][1] + n[2][2][1]) > 0))
      {
	// it matches T1
	return true;
      }
  
    // T2
    if(
       ((n[1][1][1] == OBJECT) && (n[1][2][1] == OBJECT))
       &&
       (n[0][0][0] + n[1][0][0] + n[2][0][0] +
	n[0][0][1] + n[1][0][1] + n[2][0][1] +
	n[0][0][2] + n[1][0][2] + n[2][0][2] == 0)
       && 
       (n[0][1][0] + n[1][1][0] + n[2][1][0] + n[0][1][1] + n[2][1][1] + n[0][1][2] + n[1][1][2] + n[2][1][2] + 
	n[0][2][0] + n[1][2][0] + n[2][2][0] + n[0][2][1] + n[2][2][1] + n[0][2][2] + n[1][2][2] + n[2][2][2] > 0))
      {
	// it matches T2
	return true;
      }

    // T3
    if(
       ((n[1][1][1] == OBJECT) && (n[1][2][0] == OBJECT)) 
       &&
       ((n[0][0][0] + n[1][0][0] + n[2][0][0] + 
	 n[0][0][1] + n[1][0][1] + n[2][0][1] + 
	 n[0][0][2] + n[1][0][2] + n[2][0][2] + 
	 n[0][1][2] + n[1][1][2] + n[2][1][2] + 
	 n[0][2][2] + n[1][2][2] + n[2][2][2]) == 0)
       &&
       ((n[0][1][0] + n[0][2][0] + n[2][1][0] + n[2][2][0] +
	 n[0][1][1] + n[0][2][1] + n[2][1][1] + n[2][2][1]) > 0))
      {
	// it matches T3
	return true;
      }
  
    // T4
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) && 
	(n[1][2][1] == OBJECT))
       &&
       ((n[1][0][1] + n[0][0][2] + n[1][0][2] + n[2][0][2] + n[1][1][2]) == 0)
       &&
       ((n[0][0][1] * n[0][1][2]) == 0)
       &&
       ((n[2][0][1] * n[2][1][2]) == 0))
      {
	// it matches T4
	return true;
      }
  
    // T5
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) && 
	(n[1][2][1] == OBJECT) && (n[2][0][2] == OBJECT))
       &&
       ((n[1][0][1] + n[0][0][2] + n[1][0][2] + n[1][1][2]) == 0)
       && 
       ((n[0][0][1] * n[0][1][2]) == 0)
       &&
       ((n[2][0][1] + n[2][1][2]) == OBJECT))
      {
	// matches T5
	return true;
      }

    // T6
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) && 
	(n[1][2][1] == OBJECT) && (n[0][0][2] == OBJECT))
       &&
       ((n[1][0][1] + n[1][0][2] + n[2][0][2] + n[1][1][2]) == 0)
       &&
       ((n[0][0][1] + n[0][1][2]) == OBJECT)
       &&
       ((n[2][0][1] * n[2][1][2]) == 0))
      {
	// matches T6
	return true;
      }

    // T7
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) &&
	(n[1][2][1] == OBJECT) && (n[2][1][1] == OBJECT))
       &&
       ((n[1][0][1] + n[0][0][2] + n[1][0][2] + n[1][1][2]) == 0)
       &&
       ((n[0][0][1] * n[0][1][2]) == 0))
      {
	return true;
      }

    // T8
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) && 
	(n[1][2][1] == OBJECT) && (n[0][1][1] == OBJECT))
       &&
       ((n[1][0][1] + n[1][0][2] + n[2][0][2] + n[1][1][2]) == 0)
       &&
       ((n[2][0][1] * n[2][1][2]) == 0))
      {
	return true;
      }

    // T9
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) &&
	(n[1][2][1] == OBJECT) && (n[2][1][1] == OBJECT) &&
	(n[0][0][2] == OBJECT))
       &&
       ((n[1][0][1] + n[1][0][2] + n[1][1][2]) == 0) 
       &&
       ((n[0][0][1] + n[0][1][2]) == OBJECT))
      {
	return true;
      }

    // T10
    if(
       ((n[1][1][1] == OBJECT) && (n[1][1][0] == OBJECT) &&
	(n[0][1][1] == OBJECT) && (n[1][2][1] == OBJECT) &&
	(n[2][0][2] == OBJECT))
       &&
       ((n[1][0][1] + n[1][0][2] + n[1][1][2]) == 0) 
       &&
       ((n[2][0][1] + n[2][1][2]) == OBJECT))
      {
	return true;
      }

    // T11
    if(
       ((n[1][1][1] == OBJECT) && (n[2][1][0] == OBJECT) &&
	(n[1][2][0] == OBJECT))
       &&
       ((n[0][0][0] + n[1][0][0] + 
	 n[0][0][1] + n[1][0][1] + 
	 n[0][0][2] + n[1][0][2] + n[2][0][2] +
	 n[0][1][2] + n[1][1][2] + n[2][1][2] +
	 n[0][2][2] + n[1][2][2] + n[2][2][2]) == 0))
      {
	return true;
      }

    // T12
    if(
       ((n[1][1][1] == OBJECT) && (n[0][1][0] == OBJECT) && 
	(n[1][2][0] == OBJECT))
       &&
       ((n[1][0][0] + n[2][0][0] + 
	 n[1][0][1] + n[2][0][1] + 
	 n[0][0][2] + n[1][0][2] + n[2][0][2] +
	 n[0][1][2] + n[1][1][2] + n[2][1][2] +
	 n[0][2][2] + n[1][2][2] + n[2][2][2]) == 0))
      {
	return true;
      }

    // T13
    if(
       ((n[1][1][1] == OBJECT) && (n[1][2][0] == OBJECT) && 
	(n[2][2][1] == OBJECT))
       &&
       ((n[0][0][0] + n[1][0][0] + n[2][0][0] + 
	 n[0][0][1] + n[1][0][1] + n[2][0][1] + 
	 n[0][0][2] + n[1][0][2] + n[2][0][2] + 
	 n[0][1][2] + n[1][1][2] + 
	 n[0][2][2] + n[1][2][2]) == 0))
      {
	return true;
      }

    // T14
    if(
       ((n[1][1][1] == OBJECT) && (n[1][2][0] == OBJECT) &&
	(n[0][2][1] == OBJECT))
       &&
       ((n[0][0][0] + n[1][0][0] + n[2][0][0] + 
	 n[0][0][1] + n[1][0][1] + n[2][0][1] + 
	 n[0][0][2] + n[1][0][2] + n[2][0][2] + 
	 n[1][1][2] + n[2][1][2] + 
	 n[1][2][2] + n[2][2][2]) == 0))
      {
	return true;
      }

    return false;
  }



  // transform neighborhood from a different direction
  bool VThinning::TransformNeighborhood(unsigned char n[3][3][3], 
					char direction, 
					unsigned char USn[3][3][3]) 
  {
    char i, j, k;
    unsigned char tmp[3][3][3];
  
    switch(direction) {
    case 0:  //UP_SOUTH = 0, 
      // just copy
      memcpy(USn, n, 27*sizeof(unsigned char));
      break;
    case 1:  //NORT_EAST,
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    tmp[j][2-i][k] = n[i][j][k];
	  }
	}
      }
      // 2
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[2-k][j][i] = tmp[i][j][k];
	  }
	}
      }
      break;
    case 2: // DOWN_WEST, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    tmp[2-j][i][k] = n[i][j][k];
	  }
	}
      }
      // 2
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[i][j][2-k] = tmp[i][j][k];
	  }
	}
      }
      break;
    case 3: // SOUTH_EAST,
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[2-k][j][i] = n[i][j][k];
	  }
	}
      }
      break;
    case 4: //UP_WEST, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[2-j][i][k] = n[i][j][k];
	  }
	}
      }
      break;
    case 5: // DOWN_NORTH, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    tmp[i][2-j][k] = n[i][j][k];
	  }
	}
      }
      // 2
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[i][j][2-k] = tmp[i][j][k];
	  }
	}
      }
      break;
    case 6: //SOUTH_WEST,
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[k][j][2-i] = n[i][j][k];
	  }
	}
      }
      break;
    case 7: //UP_NORTH, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[i][2-j][k] = n[i][j][k];
	  }
	}
      }
      break;
    case 8: // DOWN_EAST, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    tmp[j][2-i][k] = n[i][j][k];
	  }
	}
      }
      // 2
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[i][j][2-k] = tmp[i][j][k];
	  }
	}
      }
      break;
    case 9: // NORT_WEST,
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    tmp[2-j][i][k] = n[i][j][k];
	  }
	}
      }
      // 2
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[k][j][2-i] = tmp[i][j][k];
	  }
	}
      }
      break;
    case 10: // UP_EAST, 
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[j][2-i][k] = n[i][j][k];
	  }
	}
      }
      break;
    case 11: //  DOWN_SOUTH,
      // 1
      for(k=0; k < 3; k++) {
	for(j=0; j < 3; j++) {
	  for(i=0; i < 3; i++) {
	    USn[i][j][2-k] = n[i][j][k];
	  }
	}
      }
      break;
    }
  
    return true; 
  }

  bool VThinning::markBoundaryInDirection(unsigned char *vol, int L, int M, int N,
					  char direction)
  {
  
    int slsz = L*M;
    int idx;
    int i, j, k;
  
    // neighbor index in 18 directions (only first 6 used)
    int nb[18] = {
      +slsz - L,  // UP_SOUTH,   0
      +L + 1,   // NORT_EAST     1
      -slsz - 1,  // DOWN_WEST,  2 

      -L + 1,   // SOUTH_EAST    3
      +slsz - 1,  // UP_WEST,    4
      -slsz + L,  // DOWN_NORTH, 5 

      -L - 1,   // SOUTH_WEST    6
      +slsz + L,  // UP_NORTH,   7
      -slsz + 1,  // DOWN_EAST,  8

      +L - 1,   // NORT_WEST     9
      +slsz + 1,  // UP_EAST,   10
      -slsz - L  // DOWN_SOUTH  11

      +slsz, // UP              12
      -slsz, // DOWN            13
      +1,  // EAST,             14
      -1,  // WEST,             15
      +L,  // NORTH,            16
      -L,  // SOUTH,            17
    
    };
  
    for(k=1; k < (N-1); k++) {
      for(j=1; j < (M-1); j++) {
	for(i=1; i < (L-1); i++) {
	  idx = k*slsz + j*L + i;
	
	  if((vol[idx] == OBJECT) && (vol[idx + nb[direction]] == 0)) {
	    vol[idx] = D_BORDER;
	  }
	}
      }
    }
  
    return true;
  }

  int VThinning::Execute(unsigned char *vol,int SizeX,int  SizeY,int  SizeZ,vector<VPoint> *pExtractedIndexes=NULL)
  {
    //unsigned char *vol;
    int sz, slsz, idx;
    int nrDel;
    char dir;
    unsigned char nb[3][3][3];
    unsigned char USn[3][3][3];
  
    int L,M,N;
    int nrPasses;
    int i, j, k;
    int  maxnsp;
    bool canBeDeleted = false;

    L=SizeX;
    M=SizeY;
    N=SizeZ;  
    sz = L*M*N;
    slsz = L*M;

    // initialize global neighbors array
    // lower plane
    volNeighbors[0] = (-slsz -L -1);
    volNeighbors[1] = (-slsz -L +0);
    volNeighbors[2] = (-slsz -L +1);
    volNeighbors[3] = (-slsz +0 -1);
    volNeighbors[4] = (-slsz +0 +0);
    volNeighbors[5] = (-slsz +0 +1);
    volNeighbors[6] = (-slsz +L -1);
    volNeighbors[7] = (-slsz +L +0);
    volNeighbors[8] = (-slsz +L +1);
    // same plane
    volNeighbors[9]  = (+0 -L -1);
    volNeighbors[10] = (+0 -L +0);
    volNeighbors[11] = (+0 -L +1);
    volNeighbors[12] = (+0 +0 -1);
    volNeighbors[13] = (+0 +0 +0);
    volNeighbors[14] = (+0 +0 +1);
    volNeighbors[15] = (+0 +L -1);
    volNeighbors[16] = (+0 +L +0);
    volNeighbors[17] = (+0 +L +1);
    // upper plane
    volNeighbors[18] = (+slsz -L -1);
    volNeighbors[19] = (+slsz -L +0);
    volNeighbors[20] = (+slsz -L +1);
    volNeighbors[21] = (+slsz +0 -1);
    volNeighbors[22] = (+slsz +0 +0);
    volNeighbors[23] = (+slsz +0 +1);
    volNeighbors[24] = (+slsz +L -1);
    volNeighbors[25] = (+slsz +L +0);
    volNeighbors[26] = (+slsz +L +1);



    // set all object voxels to OBJECT
    int nobjects = 0;
    for(idx=0; idx < sz; idx++) {
      if(vol[idx] != 0) { vol[idx] = OBJECT; nobjects++; }
    }

    // std::cout << "nobjects: " << nobjects << std::endl;

    nrDel = 1;
    nrPasses = 1;
    maxnsp = 0;
  
    while(nrDel > 0) {
      // std::cout << "Pass %d...\n"<<  (int)nrPasses << std::endl;
      nrDel = 0;
      for(dir = 0; dir < 12; dir++) {
	// std::cout << "\tDir %d..." <<  (int)dir << std::endl;
	fflush(stdout);
      
	//printf("mark boundary ...\n");
	switch(dir) {
	case 0: // UP_SOUTH = 0, 
	  // UP
	  markBoundaryInDirection(vol, L, M, N, 12);
	  // SOUTH
	  markBoundaryInDirection(vol, L, M, N, 17);
	  break;
	case 1: // NORT_EAST,
	  // NOTH
	  markBoundaryInDirection(vol, L, M, N, 16);
	  // EAST
	  markBoundaryInDirection(vol, L, M, N, 14);
	  break;
	case 2: // DOWN_WEST, 
	  // DOWN
	  markBoundaryInDirection(vol, L, M, N, 13);
	  // WEST
	  markBoundaryInDirection(vol, L, M, N, 15);
	  break;
	case 3: //  SOUTH_EAST,
	  // SOUTH
	  markBoundaryInDirection(vol, L, M, N, 17);
	  // EAST
	  markBoundaryInDirection(vol, L, M, N, 14);
	  break;
	case 4: // UP_WEST, 
	  // UP
	  markBoundaryInDirection(vol, L, M, N, 12);
	  // WEST
	  markBoundaryInDirection(vol, L, M, N, 15);
	  break;
	case 5: // DOWN_NORTH, 
	  // DOWN
	  markBoundaryInDirection(vol, L, M, N, 13);
	  // NORTH
	  markBoundaryInDirection(vol, L, M, N, 16);
	  break;
	case 6: // SOUTH_WEST,
	  // SOUTH
	  markBoundaryInDirection(vol, L, M, N, 17);
	  // WEST
	  markBoundaryInDirection(vol, L, M, N, 15);
	  break;
	case 7: // UP_NORTH, 
	  // UP
	  markBoundaryInDirection(vol, L, M, N, 12);
	  // NORTH
	  markBoundaryInDirection(vol, L, M, N, 16);
	  break;
	case 8: // DOWN_EAST, 
	  // DOWN
	  markBoundaryInDirection(vol, L, M, N, 13);
	  // EAST
	  markBoundaryInDirection(vol, L, M, N, 14);
	  break;
	case 9: //  NORT_WEST,
	  // NORTH
	  markBoundaryInDirection(vol, L, M, N, 16);
	  // WEST
	  markBoundaryInDirection(vol, L, M, N, 15);
	  break;
	case 10: // UP_EAST, 
	  // UP
	  markBoundaryInDirection(vol, L, M, N, 12);
	  // EAST
	  markBoundaryInDirection(vol, L, M, N, 14);
	  break;
	case 11: // DOWN_SOUTH,
	  // DOWN
	  markBoundaryInDirection(vol, L, M, N, 13);
	  // SOUTH
	  markBoundaryInDirection(vol, L, M, N, 17);
	  break;
	}
	//printf("checking each border voxel ...\n");
	// check each boundary point and remove it if itmacthes a template
	for(k=1; k < (N-1); k++) {
	  for(j=1; j < (M-1); j++) {
	    for(i=1; i < (L-1); i++) {
	      idx = k*slsz + j*L + i;
	    
	      if(vol[idx] == D_BORDER) {
		// copy neighborhood into buffer
		//printf("copy neighborhood...\n");
		CopyNeighborhoodInBuffer(vol, L, M, N, idx, nb);
	      
		TransformNeighborhood(nb, dir, USn);
		//printf("check...\n");
		if(MatchesATemplate(USn)) {		  
		  // delete the point
		  // can be removed
		  vol[idx] = SIMPLE;
		  nrDel++;
		}
	      }
	    }
	  }
	}
      
	// reset all object voxels to OBJECT
	for(idx=0; idx < sz; idx++) {
	  // delete simple points
	  if(vol[idx] == SIMPLE) vol[idx] = 0;
	  if(vol[idx] != 0) vol[idx] = OBJECT;
	}
	// std::cout << "done.\n" << std::endl;
      }
      // std::cout << "Number of deleted voxels in pass %d: %d.\n" << (int)nrPasses << " " <<  (int)nrDel << std::endl;
      nrPasses++;
    }

    if( pExtractedIndexes )
      pExtractedIndexes->clear();
    VPoint p;
    for(k=1; k < (N); k++) {
      for(j=1; j < (M); j++) {
	for(i=1; i < (L); i++) {
	  idx = k*slsz + j*L + i;
	  if(idx>L*M*(N-1))vol[idx] =0;
	  if(vol[idx] ==OBJECT)
	    {
	      p.m_comp[0] = i;
	      p.m_comp[1] = j;
	      p.m_comp[2] = k;
	      if(pExtractedIndexes )
		pExtractedIndexes->push_back(p);
	    }
	}
      }
    }
    return 0;
  }

}
#endif
