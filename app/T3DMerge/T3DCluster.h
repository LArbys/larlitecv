#ifndef __T3DCluster_H__
#define __T3DCluster_H__

/* =================================================================
 * T3DCluster: track-3d cluster
 * 
 * This is the representation of a track that we use to do our
 * 3D cluster routines on.
 *
 * We define a track to be an ordered list of 3D points.
 * We will be doing a lot of nearest line, nearest point calculations
 * using a kDtree. To save time, we also define a bounding box in order
 * to do broad, first pass collision/overlap testing.
 * 
 * we will make use of a lot of the larlite::geotools and
 * larcv::ANN code base.
 * =================================================================*/

#include <vector>

// larlitecv/BasicTool
#include "GeoAlgo/GeoAABox.h"

#include "DataFormat/Pixel2DCluster.h"

namespace larlitecv {

  typedef std::vector<double> Point_t;
  
  class T3DCluster {

  public:
    class Builder; // we use the builder pattern to make sure this class is built ok
    T3DCluster() {};
  protected:    

    std::vector< Point_t > m_path;
    std::vector< std::vector<double> > m_dir;
    geoalgo::AABox m_bbox;
    float m_ave_stepsize;
    
  public:
    T3DCluster( const std::vector<Point_t>& path );
    T3DCluster( const std::vector<Point_t>& path, const std::vector< std::vector<double> >& pathdir, const geoalgo::AABox& bbox );
    virtual ~T3DCluster() {};

    bool overlaps( const T3DCluster& rhs ) const;
    float getAveStepSize() const { return m_ave_stepsize; };
    const std::vector< Point_t >& getPath() const { return m_path; };
    const std::vector< std::vector<double>  >& getPathDir() const { return m_dir; };
    const geoalgo::AABox& getBBox() const { return m_bbox; };
    void reverse();
    void append( const T3DCluster& end );
    void makePathDir();
    void updateBBox();    
    int pathSize() { return m_path.size(); };

    std::vector<larcv::Pixel2DCluster> getPixelsFromImages( const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							    const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							    const float stepsize );
						 
  };

  class T3DCluster::Builder {

  protected:
    
    std::vector<Point_t> path;
    
  public:
    
    Builder() {};

    Builder& setPath( const std::vector<Point_t>& path );
    Builder& addPoint( const Point_t& pt );
    int pathSize() { return path.size(); };

    T3DCluster build();

    void clear();

    
  };

}

#endif
