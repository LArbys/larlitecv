#ifndef __SEARCHALGOBASE_H__
#define __SEARCHALGOBASE_H__

// lcv
#include "Base/PSet.h"

// llcv
#include "LLCVBase/llcv_base.h"
#include "DetachedCandidate.h"

// locv
#include "LArOpenCV/Core/ImageMeta.h"
#include "LArOpenCV/ImageCluster/AlgoData/Vertex.h"
#include "LArOpenCV/ImageCluster/Base/ImageClusterFMWKInterface.h"
#include "LArOpenCV/ImageCluster/AlgoClass/LArPlaneGeo.h"

namespace llcv {

  
  //
  // Base class for searching for detached stuff in the shower image
  // user should inherit from me and fill out Configure & _Search_ methods
  // to return a vector of DetachedCluster
  //

  class SearchAlgoBase : public llcv_base {

  public:

  SearchAlgoBase(std::string name="SearchAlgoBase")
    : llcv_base(name)
      ,_name(name) 
    {}
    
    virtual ~SearchAlgoBase() {}

    const std::string& Name() const;

    virtual void Configure(const larcv::PSet &pset) = 0;

    std::vector<llcv::DetachedCandidate>
      Search(larocv::data::Vertex3D& vtx3d,
	     const std::vector<std::tuple<cv::Mat,larocv::ImageMeta> >& adc_mat_meta_v,
	     const std::vector<std::tuple<cv::Mat,larocv::ImageMeta> >& shr_mat_meta_v);
    
  protected:
    
    virtual std::vector<llcv::DetachedCandidate>
      _Search_(const larocv::data::Vertex3D& vtx3d,
	       std::vector<cv::Mat>& adc_mat_v,
	       std::vector<cv::Mat>& shr_mat_v,
	       const std::vector<larocv::ImageMeta>& meta_v) = 0;
    
    std::string _name;

  private:
    larocv::LArPlaneGeo _geo;
      
  };
}
#endif
