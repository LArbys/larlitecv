#ifndef __DBSCAN_H__
#define __DBSCAN_H__

#include <vector>
#include <cmath>
#include <cstddef>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

#include "Point.h"

namespace llcv {

  class DBSCAN {
    
  public:    
  
  DBSCAN() : m_pointSize(0), m_minPoints(0), m_epsilon(-1), _configured(false) {}
    
    ~DBSCAN(){}
    
    int run();
    
    std::vector<int> calculateCluster (const Point& point) const;
    int expandCluster(Point& point, int clusterID);

    inline float calculateDistance(const Point& pointCore, const Point& pointTarget) const;
    inline bool  calculateEquals(const Point& pt1, const Point& pt2) const;
    
    int getTotalPointSize()     const { return m_pointSize; }
    int getMinimumClusterSize() const { return m_minPoints; }
    int getEpsilonSize()        const { return m_epsilon;   }
    const std::vector<Point>& Points() const { return m_points; }

    void Configure(unsigned int minPts, float eps);
    void Reset(const std::vector<Point>& points);

  private:
    
    
    std::vector<Point> m_points;
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;
    bool _configured;

  };
  
}
#endif
