#ifndef __BOUNDARY_END_PT__
#define __BOUNDARY_END_PT__

namespace larlitecv {

  class BoundaryEndPt {
  public:

    typedef enum { kUndefined=-1, kTop, kBottom, kUpstream, kDownstream, kAnode, kCathode, kImageEnd, kNumEndTypes } BoundaryEnd_t;

    BoundaryEndPt() { t=0; w=0; type=kUndefined; };
    BoundaryEndPt( int t_, int w_, BoundaryEnd_t type_=kUndefined ) { t=t_; w=w_; type=type_; };
    virtual ~BoundaryEndPt() {}; 
    void operator= ( const BoundaryEndPt& rhs ) {
      t = rhs.t;
      w = rhs.w;
      type = rhs.type;
    };
    
    int t;
    int w;
    BoundaryEnd_t type;
    
  };
}

#endif
