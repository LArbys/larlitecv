#ifndef __STOPMU_ALGO_TYPES__
#define __STOPMU_ALGO_TYPES__

#include <array>
#include <vector>
#include <list>

namespace larlitecv {
  
  typedef std::array<float,2>  Vec2D_t;
  typedef std::array<float,3>  Vec3D_t;
  typedef std::vector<float>   Point3D_t;
  typedef std::vector<int>     Point2D_t;
  typedef std::vector< Vec2D_t > PlaneVec2D_t;
  typedef std::list< PlaneVec2D_t > Vec2DList_t;
  typedef std::list< Vec3D_t > Vec3DList_t;
  
}

#endif
