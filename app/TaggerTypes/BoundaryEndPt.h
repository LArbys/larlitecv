#ifndef __BOUNDARY_END_PT__
#define __BOUNDARY_END_PT__

#include "BoundaryMuonTaggerTypes.h"
#include <vector>
#include <stdexcept>
#include <sstream>

namespace larlitecv {

  class BoundaryEndPt {
  public:

    
    BoundaryEndPt() { row=0; col=0; type=larlitecv::kUndefined; dir.resize(2,0.0); };
    
    BoundaryEndPt( int row_, int col_, BoundaryEnd_t type_=larlitecv::kUndefined ) { 
      row=row_; 
      col=col_; 
      type=type_; 
      dir.resize(2,0);
      if ( row_<0 || col_<0 ) {
	std::stringstream msg;
	msg << __FILE__ << ":" << __LINE__ << " invalues row (" << row << ") and/or col (" << col << ") values provided" << std::endl;
	throw std::runtime_error(msg.str());
      }
    };

    virtual ~BoundaryEndPt() {}; 

    void operator= ( const BoundaryEndPt& rhs ) {
      row = rhs.row;
      col = rhs.col;
      type = rhs.type;
      dir = rhs.dir;
    };
    
    int row;
    int col;
    BoundaryEnd_t type;
    std::vector<float> dir;

    int getrow() const { return row; }
    int getcol() const { return col; }
    
  };
}

#endif
