#ifndef __VPOINT_H__
#define __VPOINT_H__

namespace llcv {
  class VPoint {
  public:
    VPoint() {}
    VPoint(const VPoint &p){
      for( int i = 0; i < 3; ++i ){
	m_comp[i] = p.m_comp[i];
      }
    }

    ~VPoint() {}
    
    int operator[](const int i) { return m_comp[i]; }

    bool operator==(const VPoint &rhs) const {
      for( int i = 0; i < 3; ++i ) {
	if( m_comp[i] != rhs.m_comp[i] ) 
	  { return false; }
      }
      return true;
    }

  public:
    int m_comp[3];

  };


  const int offset6[26][3] = { 
    {0, -1, 0}, 
    {-1, 0, 0}, 
    {1, 0, 0},
    {0, 1, 0}, 
    {0, 0, -1}, 
    {0, 0, 1}, 
  };

  const int offset26[26][3] = { 
    {-1, -1, 0}, 
    {0, -1, 0}, 
    {1, -1, 0},
    {-1, 0, 0}, 
    {1, 0, 0},
    {-1, 1, 0}, 
    {0, 1, 0}, 
    {1, 1, 0},

    {-1, -1, -1}, 
    {0, -1, -1}, 
    {1, -1, -1},
    {-1, 0, -1}, 
    {0, 0, -1}, 
    {1, 0, -1},
    {-1, 1, -1}, 
    {0, 1, -1}, 
    {1, 1, -1},

    {-1, -1, 1}, 
    {0, -1, 1}, 
    {1, -1, 1},
    {-1, 0, 1}, 
    {0, 0, 1}, 
    {1, 0, 1},
    {-1, 1, 1}, 
    {0, 1, 1}, 
    {1, 1, 1}
  };
}
#endif
