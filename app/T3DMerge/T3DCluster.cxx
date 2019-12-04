#include "T3DCluster.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TPrincipal.h"

#include "TaggerTypes/Path2Pixels.h"

namespace larlitecv {

  // =============================================================================
  // T3DCluster class
  
  T3DCluster::T3DCluster( const std::vector<T3DCluster::Point_t>& p, const std::vector< std::vector<double> >& pathdir, const ::larlite::geoalgo::AABox& bbox )
    : m_path(p), m_dir(pathdir), m_bbox(bbox), m_parent(NULL), m_is_primary(false)  {
    updatePCA();
  }

  T3DCluster::T3DCluster( const std::vector<T3DCluster::Point_t>& p )
    : m_path(p), m_parent(NULL), m_is_primary(false)  {
    update();
  }

  void T3DCluster::update() {
    makePathDir();
    updateBBox();
    updatePCA();    
  }
    
  bool T3DCluster::overlaps( const T3DCluster& rhs ) const {
    // we test to see if any of the 8 points from the rhs box is inside out box and vice versa
    bool overlaps = false;
    // // rhs corners
    // for (int i=0; i<2; i++) {
    //   for (int j=0; j<2; j++) {
    // 	for (int k=0; k<2; k++) {
    // 	  T3DCluster::Point_t corner;
    // 	  corner.resize(3);
    // 	  corner[0] = (i==0) ? rhs.m_bbox.Min()[0] : rhs.m_bbox.Max()[0];
    // 	  corner[1] = (j==0) ? rhs.m_bbox.Min()[1] : rhs.m_bbox.Max()[1];
    // 	  corner[2] = (k==0) ? rhs.m_bbox.Min()[2] : rhs.m_bbox.Max()[2];

    // 	  overlaps = true;
    // 	  for (int v=0; v<3; v++) {
    // 	    if ( corner[v] < m_bbox.Min()[v] || corner[v] > m_bbox.Max()[v] )
    // 	      overlaps = false;
    // 	  }
    // 	  //overlaps = m_bbox.Contain( corner );
    // 	  if ( overlaps )
    // 	    break;
    // 	}
    // 	if ( overlaps )
    // 	  break;
    //   }
    //   if ( overlaps )
    // 	break;
    // }
    
    // for (int i=0; i<2; i++) {
    //   for (int j=0; j<2; j++) {
    // 	for (int k=0; k<2; k++) {
    // 	  T3DCluster::Point_t corner;
    // 	  corner.resize(3);
    // 	  corner[0] = (i==0) ? m_bbox.Min()[0] : m_bbox.Max()[0];
    // 	  corner[1] = (j==0) ? m_bbox.Min()[1] : m_bbox.Max()[1];
    // 	  corner[2] = (k==0) ? m_bbox.Min()[2] : m_bbox.Max()[2];
    // 	  overlaps = true;
    // 	  for (int v=0; v<3; v++) {
    // 	    if ( corner[v] < rhs.m_bbox.Min()[v] || corner[v] > rhs.m_bbox.Max()[v] )
    // 	      overlaps = false;
    // 	  }
    // 	  //overlaps = rhs.m_bbox.Contain( corner );
    // 	  if ( overlaps )
    // 	    break;
    // 	}
    // 	if ( overlaps )
    // 	  break;
    //   }
    //   if ( overlaps )
    // 	break;
    // }

    return ( (m_bbox.Min()[0] <= rhs.m_bbox.Max()[0] && m_bbox.Max()[0] >= rhs.m_bbox.Min()[0]) &&
	     (m_bbox.Min()[1] <= rhs.m_bbox.Max()[1] && m_bbox.Max()[1] >= rhs.m_bbox.Min()[1]) &&
	     (m_bbox.Min()[2] <= rhs.m_bbox.Max()[2] && m_bbox.Max()[2] >= rhs.m_bbox.Min()[2]) );
    
    //return overlaps;
  }
  
  void T3DCluster::reverse() {
    // we reverse the sequence of path and dir list. BBox should not be affected.
    std::reverse( m_path.begin(), m_path.end() );
    update();
  }

  void T3DCluster::append( const T3DCluster& end ) {
    for ( auto const& pt : end.getPath() ) {
      m_path.push_back( pt );
    }
    update();
  }
  
  void T3DCluster::makePathDir() {
    if ( m_path.size()<2 )
      return;
    m_dir.clear();
    m_ave_stepsize = 0;
    int nsteps = 0;
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
      m_dir.push_back( segdir );
      m_ave_stepsize += dist;
      nsteps++;
    }
    m_ave_stepsize /= float(nsteps);
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

  void T3DCluster::updatePCA() {
    m_pc_v.clear();
    m_pc_bounds.clear();
    
    TPrincipal pca(3,"ND");
    Double_t x[3];
    for ( auto const& pt : getPath() ) {
      for (int i=0; i<3; i++) x[i] = pt[i];
      pca.AddRow( x );
    }
    
    pca.MakePrincipals();

    m_mean.resize(3,0);
    for ( int i=0; i<3; i++) {
      m_mean[i] = (*(pca.GetMeanValues()))[i];
    }

    for (int ipca=0; ipca<3; ipca++) {
      TVectorD eigv(3);
      eigv = TMatrixDColumn_const( *(pca.GetEigenVectors()), ipca );
      std::vector<double> vec(3);
      for (int v=0; v<3; v++) {
	vec[v] = eigv[v];
      }
      m_pc_v.push_back(vec);
      
      // orient

      // // bounds
      std::vector<double> bounds(2,0);
      const std::vector<double>& start = getPath().front();
      const std::vector<double>& end   = getPath().back();
      bounds[0] = getPCvalue( start, ipca );
      bounds[1] = getPCvalue( end, ipca );      
      
      if ( bounds[1]<0 ) {
      	for (int v=0; v<3; v++)
      	  m_pc_v[ipca][v] *= -1.0; // reorient primary PCA vector to go from start to end direction
      	bounds[0] = getPCvalue( start, 0 );
      	bounds[1] = getPCvalue( end, 0 );
      }
      m_pc_bounds.push_back(bounds);
    }
  }

  void T3DCluster::addDaughter( T3DCluster* daughter ) {
    daughter->m_parent = this;
    daughter->m_is_primary = false;
    m_daughters.push_back( daughter );
  }

  double T3DCluster::getPCvalue( const std::vector<double>& pt, int ipca) const {
    double cosdir = 0;
    std::vector<double> dx(3);
    for (int i=0; i<3; i++) {
      dx[i] = pt[i]-m_mean[i];
      cosdir += dx[i]*m_pc_v.at(ipca)[i];
    }
    return cosdir;
  }
  // =============================================================================
  // T3DCluster Builder

  void T3DCluster::Builder::clear() {
    path.clear();
  }
  
  T3DCluster::Builder& T3DCluster::Builder::setPath( const std::vector<T3DCluster::Point_t>& p ) {
    for ( auto const& pt : p )
      addPoint( pt );
    return *this;
  }

  T3DCluster::Builder& T3DCluster::Builder::addPoint( const T3DCluster::Point_t& pt ) {
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

  T3DCluster T3DCluster::Builder::build() {
    if ( path.size()==0 ) {
      throw std::runtime_error( "Empty Path in T3DCluster!" );
    }
    T3DCluster cluster( path );
    return cluster;
  }

  

}
