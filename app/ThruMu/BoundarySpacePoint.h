#ifndef __BOUNDARYSPACEPOINT__
#define __BOUNDARYSPACEPOINT__

/** 

	This class represents a 3D space point.
	It has an estimated 3D position (x,y,z) in detector coordinates.
	It also has BoundaryEndPts, which represent the corresponding positions in the wire plane images.
	
	Wanted to stop having vector<BoundaryEndPts> being passed everywhere.
 */

#include "BoundaryMuonTaggerTypes.h"
#include "BoundaryEndPt.h"

namespace larlitecv {

  class BoundarySpacePoint : public std::vector<BoundaryEndPt> {

  public:
  	// constructors
  	BoundarySpacePoint() : boundary_type( kUndefined ) { setup(); }; // default
  	BoundarySpacePoint( BoundaryEnd_t type ) : boundary_type(type) { setup(); }; // default with type

#ifndef __CINT__
#ifndef __CLING__
  	// move constructor
  	BoundarySpacePoint( BoundaryEnd_t type, std::vector<BoundaryEndPt>&& endpts )  // type and endpt vector
  	: std::vector<BoundaryEndPt>(std::move(endpts))
  	{
  		boundary_type = type;
  		setup();
  	};
#endif
#endif

  	virtual ~BoundarySpacePoint() {};

  	BoundaryEnd_t type() { return boundary_type; }
  	std::vector<float>& pos() { return m_pos; }
  	std::vector<float>& dir() { return m_dir; }

  protected:
  	BoundaryEnd_t boundary_type;
  	std::vector<float> m_pos;
  	std::vector<float> m_dir;
  	void setup() { 
  		m_pos.resize(3,0);
  		m_dir.resize(3,0);
  	};

  };


}


#endif
