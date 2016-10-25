#ifndef __BOUNDARY_END_PT__
#define __BOUNDARY_END_PT__

namespace larlitecv {

  class BoundaryEndPt {
  public:

    typedef enum { kUndefined=-1, kTop, kBottom, kUpstream, kDownstream, kAnode, kCathode, kImageEnd, kNumEndTypes } BoundaryEnd_t;
    
    BoundaryEndPt() { t=0; w=0; type=kUndefined; dir[0] = 0; dir[1]= 0.0; vec3[0] = vec3[1] = vec3[2] = 0.0; };
    BoundaryEndPt( int t_, int w_, BoundaryEnd_t type_=kUndefined ) { t=t_; w=w_; type=type_; dir[0]=dir[1] = 0.0; vec3[0] = vec3[1] = vec3[2] = 0.0; };
    virtual ~BoundaryEndPt() {}; 
    void operator= ( const BoundaryEndPt& rhs ) {
      t = rhs.t;
      w = rhs.w;
      type = rhs.type;
      dir[0] = rhs.dir[0];
      dir[1] = rhs.dir[1];
      for (int i=0; i<3; i++) vec3[i] = rhs.vec3[i];
    };
    
    int t;
    int w;
    BoundaryEnd_t type;
    float dir[2];
    float vec3[3];
    
    
  };
}

#endif
