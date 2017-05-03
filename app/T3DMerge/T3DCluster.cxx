#include "T3DCluster.h"

#include "TaggerTypes/Path2Pixels.h"

namespace larlitecv {

  // =============================================================================
  // T3DCluster class
  
  T3DCluster::T3DCluster( const std::vector<Point_t>& p, const std::vector< std::vector<double> >& pathdir, const geoalgo::AABox& bbox )
    : m_path(p), m_dir(pathdir), m_bbox(bbox) {
  }

  void T3DCluster::reverse() {
    // we reverse the sequence of path and dir list. BBox should not be affected.
    std::reverse( m_path.begin(), m_path.end() );
    makePathDir();
    updateBBox();    
  }

  void T3DCluster::append( const T3DCluster& end ) {
    for ( auto const& pt : end.getPath() ) {
      m_path.push_back( pt );
    }
    makePathDir();
    updateBBox();
  }
  
  void T3DCluster::makePathDir() {
    m_dir.clear();
    m_dir.resize( m_path.size()-1 );
    for (int idx=0; idx<(int)m_path.size()-1; idx++) {
      const std::vector<double> here = m_path[idx];
      const std::vector<double> next = m_path[idx+1];
      std::vector<double> segdir(3,0);
      float dist = 0.;
      for (int i=0; i<3; i++) {
	segdir[i] = next[i]-here[i];
	dist += segdir[i]*segdir[i];
      }
      dist = sqrt(dist);
      for (int i=0; i<3; i++)
	segdir[i] /= dist;
      m_dir[idx] = segdir;      
    }
  }

  void T3DCluster::updateBBox() {
    if ( m_path.size()==0 )
      return;
    
    std::vector<float> minx(3,0);
    std::vector<float> maxx(3,0);    
    for (int i=0; i<3; i++) {
      minx[i] = m_path[0][i];
      maxx[i] = m_path[0][i];
    }

    for (size_t i=1; i<m_path.size(); i++) {
      for (int v=0; v<3; v++) {
	if ( minx[v]>m_path[i][v] )
	  minx[v] = m_path[i][v];
	if ( maxx[v]<m_path[i][v] )
	  maxx[v] = m_path[i][v];
      }
    }

    m_bbox.Min( minx[0], minx[1], minx[2] );
    m_bbox.Max( maxx[0], maxx[1], maxx[2] );

    return;
  }

  std::vector<larcv::Pixel2DCluster> T3DCluster::getPixelsFromImages( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
								      const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
								      const float stepsize ) {
    return getTrackPixelsFromImages( m_path, imgs, badchimgs, thresholds, neighborhood_size, stepsize );
  }

  
  // =============================================================================
  // T3DCluster Builder

  void T3DCluster::Builder::clear() {
    path.clear();
    dir.clear();
  }
  
  T3DCluster::Builder& T3DCluster::Builder::setPath( const std::vector<Point_t>& p ) {
    for ( auto const& pt : p )
      addPoint( pt );
    return *this;
  }

  T3DCluster::Builder& T3DCluster::Builder::addPoint( const Point_t& pt ) {
    if (path.size()==0)
      path.push_back( pt );
    else {
      float dist = 0.;
      for (int i=0; i<3; i++)
	dist += (path.back()[i]-pt[i])*(path.back()[i]-pt[i]);
      dist = sqrt(dist);
      if ( dist>0 )
	path.push_back( pt );
    }
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

  T3DCluster::Builder& T3DCluster::Builder::buildDirList() {
    if ( path.size()<2 )
      return *this;
    
    dir.resize( path.size()-1 );
    for (int idx=0; idx<(int)path.size()-1; idx++) {
      const std::vector<double> here = path[idx];
      const std::vector<double> next = path[idx+1];
      std::vector<double> segdir(3,0);
      float dist = 0.;
      for (int i=0; i<3; i++) {
	segdir[i] = next[i]-here[i];
	dist += segdir[i]*segdir[i];
      }
      dist = sqrt(dist);
      for (int i=0; i<3; i++)
	segdir[i] /= dist;
      dir[idx] = segdir;      
    }

    return *this;
  }

  T3DCluster T3DCluster::Builder::build() {
    if ( path.size()==0 ) {
      throw std::runtime_error( "Empty Path in T3DCluster!" );
    }
    buildDirList();
    updateBBox();
    T3DCluster cluster( path, dir, bbox );
    return cluster;
  }

  

}
