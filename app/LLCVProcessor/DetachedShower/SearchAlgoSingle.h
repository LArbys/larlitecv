#ifndef __SEARCHALGOSINGLE_H__
#define __SEARCHALGOSINGLE_H__

// llcv
#include "SearchAlgoBase.h"

// locv
#include "LArOpenCV/ImageCluster/AlgoClass/PixelScan3D.h"

namespace llcv {

  class SearchAlgoSingle : public SearchAlgoBase {
  public:
  SearchAlgoSingle() : SearchAlgoBase("SearchAlgoSingle") {}
    ~SearchAlgoSingle() {}

    void Configure(const larcv::PSet &pset);

    std::vector<llcv::DetachedCandidate>
      _Search_(const larocv::data::Vertex3D& vtx3d,
	       std::vector<cv::Mat>& adc_mat_v,
	       std::vector<cv::Mat>& shr_mat_v,
	       const std::vector<larocv::ImageMeta>& meta_v);
    
  private:
    float _shower_frac;  
    float _shower_size;  
    bool _mask_vertex;

    larocv::PixelScan3D _PixelScan3D;
    
    inline bool CompareAsses(const std::array<size_t,3> & a1, const std::array<size_t,3> & a2)
    {  return (true ? (( a1[0] == a2[0]) && (a1[1] == a2[1]) && (a1[2] == a2[2])) : false); }
    
  };

}
#endif
