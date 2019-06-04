#ifndef __INTERIMAGEMANAGER_H__
#define __INTERIMAGEMANAGER_H__

//llcv
#include "LLCVBase/llcv_err.h"
#include "LLCVBase/llcv_base.h"

//me
#include "InterToolTypes.h"
#include "InterImage.h"

//locv

//lcv
#include "DataFormat/Image2D.h"
#include "LArbysImageMaker.h"

//cpp
#include <map>

namespace llcv {

  class InterModule;
  class InterMichel;
  class InterPMT;

  class InterDriver;

  class InterImageManager : public llcv_base {
    
    friend class InterModule;
    friend class InterMichel;
    friend class InterPMT;

    friend class InterDriver;
    
  public:
  InterImageManager(std::string name="InterImageManager") : 
    llcv_base(name), 
      _name(name) ,
      _iimg_v(nullptr)
      { Reset(); }

    ~InterImageManager(){}

  public:

    template <class T>
      std::vector<T*> Image(llcv::InterImageType iitype, int cropx, int cropy);

    template <class T>
      std::vector<T*> RawImage(llcv::InterImageType iitype );
    
  private:

    std::string _name;

    void SetImage(const std::vector<larcv::Image2D>& img_v, llcv::InterImageType iitype);
    void SetVertex(float x, float y, float z);
    void SetPixel(int row, int col, size_t plane);
    void Reset();
    void Erase();

    void SetIIT(InterImageType iitype,const std::pair<int,int>& cpair);
    void SetIIT(InterImageType iitype,int cropx, int cropy);

    void InitializeOCVImage(InterImageType iitype);

    void CropImage(int cropx, int cropy,InterImageType iitype);
    void CropImage(const std::pair<int,int>& cpair, InterImageType iitype);

    std::map<std::pair<int,int>, std::vector<InterImage> > _inter_adc_m;
    std::map<std::pair<int,int>, std::vector<InterImage> > _inter_shr_m;
    std::map<std::pair<int,int>, std::vector<InterImage> > _inter_trk_m;
    std::map<std::pair<int,int>, std::vector<InterImage> > _inter_dead_m;

    larcv::LArbysImageMaker _larmkr;
    
    std::vector<InterImage>* _iimg_v;

    std::vector<std::pair<int,int> > _vtx_pixel_v;
  };


}

#endif
