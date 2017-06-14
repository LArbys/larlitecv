#ifndef __ASTAR3DPOSTPROCESSOR_H__
#define __ASTAR3DPOSTPROCESSOR_H__

#include <vector>
#include "TaggerTypes/BMTrackCluster3D.h"

/* This function tries to make up for the failures of the AStar3DFitter.
   It takes in the tracks after the initial fitting has taken place and removes one of two tracks that are overlaid over the same region. */


namespace larlitecv {

  class AStar3DPostProcessor {

  public: 
    AStar3DPostProcessor() {};
    virtual ~AStar3DPostProcessor() {};

  std::vector< BMTrackCluster3D > separateAStarTracks(std::vector < BMTrackCluster3D >& tracks_v, double maximum_distance, const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs, const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size);

  };

}

#endif
