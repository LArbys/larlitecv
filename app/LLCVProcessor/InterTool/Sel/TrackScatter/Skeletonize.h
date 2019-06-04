#ifndef __SKELETONIZE_H__
#define __SKELETONIZE_H__

#include "VThinning.h"

#include <vector>

#include "LArOpenCV/ImageCluster/AlgoData/Vertex.h"

namespace llcv {

  class Skeletonize {
  public:
  Skeletonize() : _voxel(nullptr) {}
    ~Skeletonize() {}

    void Initialize(const std::vector<larocv::data::Vertex3D>& vtx_v,
		    float voxel_size);

    std::vector<larocv::data::Vertex3D> Run();
    void Clear();
    

  private:

    VThinning _VThinning;
    void Dump();

    int _LocationToVoxel(float x, float y, float z) const;
    larocv::data::Vertex3D _VoxelToLocation(int idx) const;
    larocv::data::Vertex3D _VoxelToLocation(int x_loc, int y_loc, int z_loc) const;
    larocv::data::Vertex3D _VoxelToLocation(const VPoint& pt) const;

  private:

    unsigned char* _voxel;
    std::vector<VPoint> _thin_voxel_v;

    float _min_x;
    float _min_y;
    float _min_z;

    float _max_x;
    float _max_y;
    float _max_z;

    int _lX;
    int _lY;
    int _lZ;

    float _voxel_size;

    int _sz;

  };

}


#endif
