#include "ContourShapeMeta.h"
#include <sstream>

namespace larlitecv {

  ContourShapeMeta::ContourShapeMeta() {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << ": Default construct should not be used. Only defined for dictionary purposes." << std::endl;
    throw std::runtime_error(msg.str());
  }
  
  ContourShapeMeta::ContourShapeMeta( const std::vector<cv::Point>& contour, const larcv::ImageMeta& meta )
    : std::vector<cv::Point>(contour),
      m_meta(meta),
      m_start( cv::Point(0,0) ),
      m_end( cv::Point(0,0) )  {
    _fill_linefit_members();
    _build_bbox();
    _get_tick_range();
  }

  void ContourShapeMeta::_fill_linefit_members() {
    cv::Vec4f out_array;

    cv::fitLine( *this, out_array, cv::DIST_L2, 0, 0.01, 0.01 );
    m_dir.resize(2,0);
    
    // loop through contour points and get min and max on projected line
    if ( out_array[1]!=0 ) {
      float norm  = sqrt(out_array[0]*out_array[0] + out_array[1]*out_array[1]);
      m_dir[0] = out_array[0]/norm;
      m_dir[1] = out_array[1]/norm;
    }
    else {
      // vertical
      float norm  = sqrt(out_array[0]*out_array[0] + out_array[1]*out_array[1]);
      m_dir[0] = out_array[0]/norm;
      m_dir[1] = 0;
    }
    
    cv::Point maxpt(0,0);
    cv::Point minpt(0,0);
    float mincos =  1e6;
    float maxcos = -1e6;
    for ( auto& pt : (*this) ) {
      float dx[2];
      dx[0] = pt.x-out_array[2];
      dx[1] = pt.y-out_array[3];
      float ptcos = 0.;
      for (int i=0; i<2; i++)
	ptcos += dx[i]*m_dir[i];
      if ( ptcos < mincos ) {
	minpt = pt;
	mincos = ptcos;
      }
      if ( ptcos > maxcos ) {
	maxpt = pt;
	maxcos = ptcos;
      }
    }
    m_start = minpt;
    m_end   = maxpt;

    // orient: start is at low y
    if ( m_start.y > m_end.y ) {
      cv::Point temp = m_start;
      m_start = m_end;
      m_end = temp;
      for (int i=0; i<2; i++)
	m_dir[i] *= -1.;
    }
    
  }
  
  void ContourShapeMeta::_build_bbox() {
    m_bbox = cv::boundingRect( *this );
  }

  void ContourShapeMeta::_get_tick_range() {
    float miny = -1;
    float maxy = -1;
    float minx = -1;
    float maxx = -1;
    for (auto& pt : *this ) {
      if ( miny<0 || pt.y<miny )
	miny = pt.y;
      if (maxy<0 || pt.y>maxy )
	maxy = pt.y;
      if ( minx<0 || pt.x<minx )
	minx = pt.x;
      if ( maxx<0 || pt.x>maxx )
	maxx = pt.x;
    }

    ybounds.resize(2,0);
    xbounds.resize(2,0);
    ybounds[0] = miny;
    ybounds[1] = maxy;
    xbounds[0] = minx;
    xbounds[1] = maxx;
  }
  
}
