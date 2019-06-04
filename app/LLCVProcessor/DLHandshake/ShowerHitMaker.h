#ifndef SHOWERHITMAKER_H
#define SHOWERHITMAKER_H

#include "LLCVBase/AnaBase.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

namespace llcv {

  class ShowerHitMaker : public AnaBase {

  public:

  ShowerHitMaker(const std::string name="ShowerHitMaker") : AnaBase(name) {}
    ~ShowerHitMaker() {}
    
    void configure(const larcv::PSet&);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();

  private:
    
    bool MakeShowerImage(const larcv::EventPixel2D* ev_pixel2d,
			 std::vector<larcv::Image2D>& shr_img_v);
    
    
  };

}

#endif
