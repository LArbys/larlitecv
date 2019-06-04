#ifndef __INTERDATAMANAGER_CXX__
#define __INTERDATAMANAGER_CXX__

#include "InterDataManager.h"

namespace llcv { 
  
  larlite::track& InterDataManager::MakeTrack() {
    _out_track_v.resize(_out_track_v.size()+1);
    return _out_track_v.back();
  }
  
  larlite::shower& InterDataManager::MakeShower() {
    _out_shower_v.resize(_out_shower_v.size()+1);
    return _out_shower_v.back();
  }
 
  larcv::PGraph& InterDataManager::MakePGraph() {
    _out_pgraph_v.resize(_out_pgraph_v.size()+1);
    return _out_pgraph_v.back();
  }

  larcv::Pixel2DCluster& InterDataManager::MakeImage(const size_t plane, const larcv::ImageMeta& meta) {
    _out_pimg_v.resize(_out_pimg_v.size()+1);
    _out_pimg_plane_v.resize(_out_pimg_plane_v.size()+1);
    _out_pimg_meta_v.resize(_out_pimg_meta_v.size()+1);

    _out_pimg_plane_v.back() = plane;
    _out_pimg_meta_v.back() = meta;

    return _out_pimg_v.back();
  }


  larcv::Pixel2DCluster& InterDataManager::MakeInteraction(const size_t plane, const larcv::ImageMeta& meta) {
    _out_pinter_v.resize(_out_pinter_v.size()+1);
    _out_pinter_plane_v.resize(_out_pinter_plane_v.size()+1);
    _out_pinter_meta_v.resize(_out_pinter_meta_v.size()+1);

    _out_pinter_plane_v.back() = plane;
    _out_pinter_meta_v.back() = meta;

    return _out_pinter_v.back();
  }

  larcv::Pixel2DCluster& InterDataManager::MakePixel2DCluster(const size_t plane, const larcv::ImageMeta& meta) {
    _out_pixel_cluster_v.resize(_out_pixel_cluster_v.size()+1);
    _out_pixel_cluster_plane_v.resize(_out_pixel_cluster_plane_v.size()+1);
    _out_pixel_cluster_meta_v.resize(_out_pixel_cluster_meta_v.size()+1);

    _out_pixel_cluster_plane_v.back() = plane;
    _out_pixel_cluster_meta_v.back() = meta;

    return _out_pixel_cluster_v.back();
  }
 
}

#endif
