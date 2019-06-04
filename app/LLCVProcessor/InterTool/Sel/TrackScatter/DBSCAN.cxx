#ifndef __DBSCAN_CXX__
#define __DBSCAN_CXX__

#include "DBSCAN.h"

#include <stdexcept>

namespace llcv {

  void DBSCAN::Configure(unsigned int minPts, float eps) {
    m_minPoints = minPts;
    m_epsilon = eps;
    _configured = true;
  }

  void DBSCAN::Reset(const std::vector<Point>& points) {
    if (!_configured) 
      throw std::runtime_error("You must call configure");
    m_points.clear();
    m_points = points;
    m_pointSize = points.size();
  }

  int DBSCAN::run() {
    int clusterID = 1;
    for(auto& point : m_points) {
      if ( point.clusterID == UNCLASSIFIED ) {
	if ( expandCluster(point, clusterID) != FAILURE ) {
	  clusterID += 1;
	}
      }
    }
    return 0;
  }

  int DBSCAN::expandCluster(Point& point, int clusterID) {

    std::vector<int> clusterSeeds = calculateCluster(point);

    if ( clusterSeeds.size() < m_minPoints ) {
      point.clusterID = NOISE;
      return FAILURE;
    }
    
    else {
      int indexCorePoint = 0;
      for(int index=0; index<(int)clusterSeeds.size(); ++index) {

	const auto seed_id = clusterSeeds[index];
	auto& seed_pt = m_points[seed_id];
	
	seed_pt.clusterID = clusterID;
	
	if (calculateEquals(seed_pt,point))
	  indexCorePoint = index;
      }
      
      clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);

      for( size_t i = 0, n = clusterSeeds.size(); i < n; ++i ) {
	const auto seed_id = clusterSeeds[i];
	auto& seed_pt = m_points[seed_id];
	auto clusterNeighbors = calculateCluster(seed_pt);
	if (clusterNeighbors.size() >= m_minPoints) {
	  for(auto neighbor_id : clusterNeighbors) {
	    auto& neighbor_pt = m_points[neighbor_id];
	    if (neighbor_pt.clusterID == UNCLASSIFIED or neighbor_pt.clusterID == NOISE) {
	      if (neighbor_pt.clusterID == UNCLASSIFIED) {
		clusterSeeds.push_back(neighbor_id);
		n = clusterSeeds.size();
	      }
	      neighbor_pt.clusterID = clusterID;
	    }
	  }
	}
      }
    } // end else
    return SUCCESS;
  }


  std::vector<int> DBSCAN::calculateCluster (const Point& point) const {

    std::vector<int> clusterIndex;
    clusterIndex.reserve((size_t) (((float)m_points.size())/2.0));
    
    for(int index=0; index<(int)m_points.size(); ++index) {
      const auto& local_point = m_points[index];
      if (calculateDistance(point, local_point) <= m_epsilon)
	clusterIndex.push_back(index);
    }
    
    return clusterIndex;
  }

  inline float DBSCAN::calculateDistance(const Point& pointCore, const Point& pointTarget ) const {
    return pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2);
  }

  inline bool DBSCAN::calculateEquals(const Point& pt1, const Point& pt2) const {
    return (pt1.x == pt2.x and pt1.y == pt2.y and pt1.z == pt2.z);
  }
  
}

#endif
