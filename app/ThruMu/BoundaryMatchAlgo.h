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
    
    BoundaryCombo() { u=0; v=0; y=0; pos.resize(3,0); };
    BoundaryCombo( unsigned short _uwire, unsigned short _vwire, unsigned short _ywire ) {
      u = _uwire;
      v = _vwire;
      y = _ywire;
      pos.resize(3,0);
    };
    virtual ~BoundaryCombo() {};
    bool operator<( const BoundaryCombo& rhs ) const {
      if ( u<rhs.u ) return true;
      if ( u==rhs.u && v<rhs.v ) return true;
      if ( u==rhs.u && v==rhs.v && y<rhs.y ) return true;
      return false;
    };
    
    unsigned short u;
    unsigned short v;
    unsigned short y;
    std::vector<float> pos;

  };

  class BoundaryMatchList {
    
  public:


    BoundaryMatchList( BoundaryMatchArrays::Boundary_t t ) { load(t); };
    virtual ~BoundaryMatchList();

    std::vector<BoundaryCombo> findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires, 
					   int urange, int vrange, int yrange, const std::vector<larcv::Image2D>& badchs, bool use_badchs );

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

    void load( BoundaryMatchArrays::Boundary_t t  );

  };

  class BoundaryMatchAlgo {

  public:
    BoundaryMatchAlgo();
    virtual ~BoundaryMatchAlgo();

    void findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires,
		     int urange, int vrange, int yrange,
		     const std::vector<larcv::Image2D>& badchs, bool use_badchs,
		     std::vector<  std::vector<BoundaryCombo>  >& boundary_combos );
    
  protected:

    BoundaryMatchList* matchlist[4];
    
  };


}

#endif

