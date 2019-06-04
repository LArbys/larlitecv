#ifndef __OBJECT3D_H__
#define __OBJECT3D_H__

#include "LArOpenCV/ImageCluster/AlgoData/Vertex.h"
#include "DataFormat/track.h"

namespace llcv {
  
  class Object3D {
  public:
    
  Object3D(const larocv::data::Vertex3D& s,
	   const std::vector<larocv::data::Vertex3D>& v) 
    : _start(s), _pts_v(v) { FillPCA(); }
    
    Object3D(const larocv::data::Vertex3D& s,
	     const std::vector<const larocv::data::Vertex3D*> & v);
    
    ~Object3D() {}

    float Width1() const { return _width1; }
    float Width2() const { return _width2; }

    float Theta()   const { return _theta;  }
    float Phi()     const { return _phi;    }
    float Opening1() const { return _oangle1; }
    float Opening2() const { return _oangle2; }

    float Length()  const { return _length; }
    float Width()   const { return ((_width1 + _width2) / 2.0);   }
    float Opening() const { return ((_oangle1 + _oangle2) / 2.0); }

    const larocv::data::Vertex3D& Start()  const { return _start;  }
    const larocv::data::Vertex3D& End()    const { return _end;    }
    const larocv::data::Vertex3D& Center() const { return _center; }
    const larocv::data::Vertex3D& Edge1()  const { return _edge1;  }
    const larocv::data::Vertex3D& Edge2()  const { return _edge2;  }

    const std::vector<float>& PCADeviation() const { return _deviation_v; }

    const larocv::data::Vertex3D& Point(size_t pid) const { return _pts_v.at(pid); }
    const std::vector<larocv::data::Vertex3D>& Points() const { return _pts_v; }

    std::vector<float> TrackDeviation(const larlite::track& trk) const;

  private:

    float _length;
    float _width1;
    float _width2;
    float _theta;
    float _phi;
    float _oangle1;
    float _oangle2;

    larocv::data::Vertex3D _start;
    larocv::data::Vertex3D _end;
    larocv::data::Vertex3D _center;

    larocv::data::Vertex3D _edge1;
    larocv::data::Vertex3D _edge2;

    std::vector<larocv::data::Vertex3D> _pts_v;
    std::vector<float> _deviation_v;

    void FillPCA();
    void FillOOBB(const std::array<float,3> mean_v, const std::array<std::array<float,3>, 3> eigen_vv);

    std::array<float,3> ToVector(const larocv::data::Vertex3D& vtx) const;
    std::array<float,3> ToVector(const TVector3& pt) const;
    
  };


  
}
#endif
