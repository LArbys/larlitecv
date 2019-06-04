#ifndef __SKELETONIZE_CXX__
#define __SKELETONIZE_CXX__

#include "Skeletonize.h"

#include <cassert>
#include <iostream>

namespace llcv {

  
  void Skeletonize::Initialize(const std::vector<larocv::data::Vertex3D>& vtx_v, float voxel_size) {
    //
    // get the absolute scale boundaries
    //
    Clear();
    
    for(size_t vid=0; vid<vtx_v.size(); ++vid) {
      const auto& vtx = vtx_v[vid];
      _min_x = std::min(_min_x,(float)vtx.x);
      _min_y = std::min(_min_y,(float)vtx.y);
      _min_z = std::min(_min_z,(float)vtx.z);

      _max_x = std::max(_max_x,(float)vtx.x);
      _max_y = std::max(_max_y,(float)vtx.y);
      _max_z = std::max(_max_z,(float)vtx.z);
    }

    float lX = (_max_x - _min_x) + 1;
    float lY = (_max_y - _min_y) + 1;
    float lZ = (_max_z - _min_z) + 1;
    
    _voxel_size = voxel_size;

    lX /= _voxel_size;
    lY /= _voxel_size;
    lZ /= _voxel_size;

    _lX = lX;
    _lY = lY;
    _lZ = lZ;

    _sz = _lX * _lY * _lZ;
    
    
    _voxel = new unsigned char[_sz];

    for(size_t xid=0; xid<(size_t)_sz; ++xid)
      _voxel[xid] = ((unsigned char)0);
    
    for(size_t vid=0; vid<vtx_v.size(); ++vid) {
      const auto& vtx = vtx_v[vid];
      auto idx = _LocationToVoxel(vtx.x,vtx.y,vtx.z);
      _voxel[idx] = ((unsigned char)1);
    }

    return;
  }
  
  int Skeletonize::_LocationToVoxel(float x, float y, float z) const {

    int x_loc = (int)(((x - _min_x) / _voxel_size));
    int y_loc = (int)(((y - _min_y) / _voxel_size));
    int z_loc = (int)(((z - _min_z) / _voxel_size));

    assert (x_loc < _lX);
    assert (y_loc < _lY);
    assert (z_loc < _lZ);

    assert (x_loc >= 0);
    assert (y_loc >= 0);
    assert (z_loc >= 0);
    
    int idx = (_lX * _lY) * z_loc + (_lX) * y_loc + x_loc;
    
    assert (idx < _sz);
    assert (idx >=  0);

    return idx;
  }

  larocv::data::Vertex3D Skeletonize::_VoxelToLocation(int idx) const {
    
    larocv::data::Vertex3D res;

    int x_loc = idx % _lX;
    int y_loc = (idx / _lX) % _lY;
    int z_loc = idx / (_lY * _lX);

    float x = x_loc*_voxel_size + _min_x;
    float y = y_loc*_voxel_size + _min_y;
    float z = z_loc*_voxel_size + _min_z;
    
    res.x = x;
    res.y = y;
    res.z = z;

    return res;
  }

  larocv::data::Vertex3D Skeletonize::_VoxelToLocation(const VPoint& pt) const {
    return _VoxelToLocation(pt.m_comp[0],pt.m_comp[1],pt.m_comp[2]);
  }


  larocv::data::Vertex3D Skeletonize::_VoxelToLocation(int x_loc, int y_loc, int z_loc) const {
    
    larocv::data::Vertex3D res;

    float x = x_loc*_voxel_size + _min_x;
    float y = y_loc*_voxel_size + _min_y;
    float z = z_loc*_voxel_size + _min_z;
    
    res.x = x;
    res.y = y;
    res.z = z;

    return res;
  }


  std::vector<larocv::data::Vertex3D> Skeletonize::Run() {
    std::vector<larocv::data::Vertex3D> res_v;

    _VThinning.Execute(_voxel,_lX,_lY,_lZ,&_thin_voxel_v);
    
    res_v.resize(_thin_voxel_v.size());

    for(size_t id = 0; id < res_v.size(); ++id) 
      res_v[id] = _VoxelToLocation(_thin_voxel_v[id]);
    
    return res_v;
  }
  
  void Skeletonize::Dump() {

    std::cout << "_min_x: " << _min_x << std::endl;
    std::cout << "_min_y: " << _min_y << std::endl;
    std::cout << "_min_z: " << _min_z << std::endl;

    std::cout << "_max_x: " << _max_x << std::endl;
    std::cout << "_max_y: " << _max_y << std::endl;
    std::cout << "_max_z: " << _max_z << std::endl;

    std::cout << "_lX: " << _lX << std::endl;
    std::cout << "_lY: " << _lY << std::endl;
    std::cout << "_lZ: " << _lZ << std::endl;

    std::cout << "_voxel_size: " << _voxel_size << std::endl;

    std::cout << "_sz: " << _sz << std::endl;

    return;
  }

  void Skeletonize::Clear() {

    if(_voxel) 
      delete[] _voxel;

    _voxel = nullptr;

    _thin_voxel_v.clear();

    _min_x = larocv::kINVALID_FLOAT;
    _min_y = larocv::kINVALID_FLOAT;
    _min_z = larocv::kINVALID_FLOAT;

    _max_x = -1.0*larocv::kINVALID_FLOAT;
    _max_y = -1.0*larocv::kINVALID_FLOAT;
    _max_z = -1.0*larocv::kINVALID_FLOAT;

    _lX = -1.0*larocv::kINVALID_INT;
    _lY = -1.0*larocv::kINVALID_INT;
    _lZ = -1.0*larocv::kINVALID_INT;

    _sz = -1.0*larocv::kINVALID_INT;
    
    _voxel_size = -1.0*larocv::kINVALID_FLOAT;

    return;
  }

}

#endif
