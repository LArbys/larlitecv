#ifndef __ContourShapeMeta__
#define __ContourShapeMeta__

#include <vector>

#include "DataFormat/ImageMeta.h"

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif


namespace larlitecv {

 class ContourShapeMeta : public std::vector<cv::Point> {
   // Wrapper around OpenCV contour
   // Stores meta data for contours
 public:

   ContourShapeMeta() {}; // shouldn't use this
   ContourShapeMeta( const std::vector<cv::Point>& contour, const larcv::ImageMeta& img );
   virtual ~ContourShapeMeta() {};

   const larcv::ImageMeta& meta() const { return m_meta; };    
   const cv::Point& getFitSegmentStart() const { return m_start; };
   const cv::Point& getFitSegmentEnd() const { return m_end; };   
   const cv::Rect&  getBBox() const  { return m_bbox; };

   std::vector<float> getEndDir() const { return m_dir; };
   std::vector<float> getStartDir() const {
     std::vector<float> reverse_dir(m_dir.size(),0);
     for (size_t i=0; i<m_dir.size(); i++) reverse_dir[i] = -1.0*m_dir[i];
     return reverse_dir;
   };
   
   
 protected:
   
   // ImageMeta
   const larcv::ImageMeta m_meta;
   
   // 2D PCA

   // Line Fit/Projected End
   std::vector<float> m_dir;
   cv::Point m_start;
   cv::Point m_end;
   void _fill_linefit_members();

   // Bounding Box (for collision detection)
   cv::Rect m_bbox;
   void _build_bbox();
   
 };
 

}

#endif
