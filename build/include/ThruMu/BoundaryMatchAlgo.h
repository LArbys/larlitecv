/**
   
   This class constains the algo to find boundary match combinations given a list of hit wires from 3 planes

 */

#ifndef __BUONDARY_MATCH_ALGO_H__
#define __BUONDARY_MATCH_ALGO_H__

#include <vector>
#include <set>
#include <map>

#include "DataFormat/Image2D.h"
#include "BoundaryMatchArrays.h"

namespace larlitecv {

  class BoundaryCombo {
  public:
    
    BoundaryCombo() { planewire.resize(3,0); pos.resize(3,0); };
    BoundaryCombo( unsigned short _uwire, unsigned short _vwire, unsigned short _ywire, float x, float y, float z ) {
      planewire.resize(3,0);
      planewire[0] = _uwire;
      planewire[1] = _vwire;
      planewire[2] = _ywire;
      pos.resize(3,0);
      pos[0] = x;
      pos[1] = y;
      pos[2] = z;
    };
    virtual ~BoundaryCombo() {};
    
    unsigned short u() const { return planewire[0]; };
    unsigned short v() const { return planewire[1]; };
    unsigned short y() const { return planewire[2]; };
    
    bool operator<( const BoundaryCombo& rhs ) const {
      if ( u()<rhs.u() ) return true;
      if ( u()==rhs.u() && v()<rhs.v() ) return true;
      if ( u()==rhs.u() && v()==rhs.v() && y()<rhs.y() ) return true;
      return false;
    };
    
    std::vector<unsigned short> planewire;
    std::vector<float> pos;

  };
  
  class BoundaryMatchList {
    
  public:


    BoundaryMatchList( BoundaryMatchArrays::Boundary_t t, BoundaryMatchArrays::MatchMode_t mode ) { load(t, mode); };
    virtual ~BoundaryMatchList();

    std::vector<BoundaryCombo> findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires, 
					   const std::vector<larcv::Image2D>& badchs, bool use_badchs );

  protected:

    typedef unsigned short idx;
    std::vector< BoundaryCombo > m_combos;
    std::map< idx, std::set<idx> > u_indices;
    std::map< idx, std::set<idx> > v_indices;
    std::map< idx, std::set<idx> > y_indices;
    idx uid_max;
    idx vid_max;
    idx yid_max;

    int* uvec;
    int* vvec;
    int* yvec;

    void load( BoundaryMatchArrays::Boundary_t t, BoundaryMatchArrays::MatchMode_t mode  );

  };

  class BoundaryMatchAlgo {

  public:
    BoundaryMatchAlgo( BoundaryMatchArrays::MatchMode_t mode );
    virtual ~BoundaryMatchAlgo();

    void findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires,
		     const std::vector<larcv::Image2D>& badchs, bool use_badchs,
		     std::vector<  std::vector<BoundaryCombo>  >& boundary_combos );
    
  protected:

    BoundaryMatchArrays::MatchMode_t fMode;
    BoundaryMatchList* matchlist[4];
    
  };


}

#endif

