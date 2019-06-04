#ifndef __POINT_H__
#define __POINT_H__
namespace llcv {
  class Point {
  public:
  Point() : clusterID(-1) {}
    ~Point() {}
    float x;
    float y;
    float z;
    int clusterID;
  };
}
#endif
