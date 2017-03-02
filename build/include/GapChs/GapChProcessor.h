#ifndef __GAPCHPROCESSOR__
#define __GAPCHPROCESSOR__

#include <vector>

// LArCV
#include "DataFormat/Image2D.h"

namespace larlitecv {

  class GapChProcessor {
  public:
    GapChProcessor() {};
    virtual ~GapChProcessor() {};

    void addEventImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );

    std::vector<larcv::Image2D> makeGapChImage();

  protected:

  	std::vector< std::vector<larcv::Image2D> > m_evimg_list;
  	std::vector< std::vector<larcv::Image2D> > m_evbadch_list;


  };

}

#endif
