#ifndef __ASTAR_NODES_2_BMTRACKCLUSTER3D_H__
#define __ASTAR_NODES_2_BMTRACKCLUSTER3D_H__

#include <vector>

#include "DataFormat/Image2D.h"
#include "DataFormat/Pixel2D.h"

#include "ThruMu/BMTrackCluster3D.h"
#include "ThruMu/AStar3DAlgo.h"


namespace larlitecv {

  BMTrackCluster3D AStarNodes2BMTrackCluster3D( const std::vector<AStar3DNode>& path, 
    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
    const std::vector<const larcv::Pixel2D*>& start_pt, const std::vector<const larcv::Pixel2D*>& end_pt,
    const int pixel_tag_neighborhood, const float link_step_size, const std::vector<float>& pixel_thresholds );


}

#endif