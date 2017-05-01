#include "T3DCluster.h"

namespace larlitecv {

  // =============================================================================
  // T3DCluster constructor
  
  T3DCluster::T3DCluster( const std::vector<Point_t>& p, const geoalgo::AABox& bbox )
    : m_path(p), m_bbox(bbox) {
  }

  // =============================================================================
  // T3DCluster Builder
  
  T3DCluster::Builder& T3DCluster::Builder::setPath( const std::vector<Point_t>& p ) {
    path = p;
    return *this;
  }

  T3DCluster::Builder& T3DCluster::Builder::addPoint( const Point_t& pt ) {
    path.push_back( pt );
    return *this;
  }

  T3DCluster::Builder& T3DCluster::Builder::updateBBox() {
    if ( path.size()==0 )
      return *this;

    std::vector<float> minx(3,0);
    std::vector<float> maxx(3,0);    
    for (int i=0; i<3; i++) {
      minx[i] = path[0][i];
      maxx[i] = path[0][i];
    }

    for (size_t i=1; i<path.size(); i++) {
      for (int v=0; v<3; v++) {
	if ( minx[v]>path[i][v] )
	  minx[v] = path[i][v];
	if ( maxx[v]<path[i][v] )
	  maxx[v] = path[i][v];
      }
    }

    bbox.Min( minx[0], minx[1], minx[2] );
    bbox.Max( maxx[0], maxx[1], maxx[2] );

    return *this;
  }

  T3DCluster T3DCluster::Builder::build() {
    if ( path.size()==0 ) {
      throw std::runtime_error( "Empty Path in T3DCluster!" );
    }
    updateBBox();
    T3DCluster cluster( path, bbox );
    return cluster;
  }

}
