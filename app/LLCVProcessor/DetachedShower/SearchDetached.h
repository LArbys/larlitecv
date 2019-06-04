#ifndef SEARCHDETACHED_H
#define SEARCHDETACHED_H

// me
#include "SearchAlgoBase.h"

//llcv
#include "LLCVBase/AnaBase.h"

//lcv
#include "DataFormat/EventPGraph.h"
#include "LArbysImageMaker.h"

//ll

namespace llcv {
  
  
  //
  // Simple Algo driver which takes pgraph and pixels, calls your algorithm
  // to search for detached pixels, if you find something, it makes a
  // pgraph and pixel2d cluster for you. Send my output to handshaker & shower reco
  //

  class SearchDetached : public AnaBase {
    
  public:
    
  SearchDetached(const std::string name="SearchDetached")
    :  AnaBase(name)
      ,_algo(nullptr) {}

    ~SearchDetached() {}
    
    void configure(const larcv::PSet& cfg);
    void initialize();
    bool process(larcv::IOManager& mgr, larlite::storage_manager& sto);
    void finalize();
    
    // larcv -> cv::Mat image convert
    const larcv::LArbysImageMaker& LArMkr() const { return _larmkr; }

    // given an image, mask out the particles
    void MaskImage(const std::vector<larcv::PGraph>& pgraph_v,
		   const std::map<larcv::PlaneID_t, std::vector<larcv::Pixel2DCluster> >& pix_m,
		   const std::map<larcv::PlaneID_t, std::vector<larcv::ImageMeta> >&   pix_meta_m,
		   std::vector<larcv::Image2D>& img_v);

    // set the search algorithm to be run -- one at a time
    void SetAlgo(SearchAlgoBase* algo) 
    { if(!algo) throw llcv_err("Invalid algo specified"); _algo = algo; }

  private:

    // given vector of detached candidates, fill the output ptrs
    void FillOutput(const std::vector<DetachedCandidate>& detached_v,
		    larcv::EventPGraph* ev_pgraph,
		    larcv::EventPixel2D* ev_pixel,
		    const std::vector<larcv::Image2D>& adc_img_v);

  private:

    SearchAlgoBase* _algo;
    
    std::string _adc_img_prod;
    std::string _shr_img_prod;
    std::string _pgraph_prod;
    std::string _pixel_prod;

    std::string _out_pgraph_prod;
    std::string _out_pixel_prod;
    
    float _search_distance;
    bool _mask_particles;

    larcv::LArbysImageMaker _larmkr;
    
  };

}

#endif
