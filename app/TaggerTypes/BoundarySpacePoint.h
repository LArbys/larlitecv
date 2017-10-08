#ifndef __BOUNDARYSPACEPOINT__
#define __BOUNDARYSPACEPOINT__

/**

	This class represents a 3D space point.
	It has an estimated 3D position (x,y,z) in detector coordinates.
	It also has BoundaryEndPts, which represent the corresponding positions in the wire plane images.

	Wanted to stop having vector<BoundaryEndPts> being passed everywhere.

	Going from boundaryendpts -> (x,y,z) needs to be better thought out. handled externally. seems dumb.
 */

#include "BoundaryMuonTaggerTypes.h"
#include "BoundaryEndPt.h"

#include <stdexcept>
#include <vector>
#include <string>

// larlite
#include "DataFormat/opflash.h"

// larcv
#include "DataFormat/ImageMeta.h"

#include "BoundaryFlashIndex.h"

namespace larlitecv {

  class BoundarySpacePoint : public std::vector<BoundaryEndPt> {

  public:
    // constructors
  BoundarySpacePoint()
    : boundary_type( kUndefined ) {
      m_pos.resize(3,0);
      m_dir.resize(3,0);
      m_empty=true;
    }; // default
  BoundarySpacePoint( BoundaryEnd_t type )
    : boundary_type(type)
    {
      m_pos.resize(3,0);
      m_dir.resize(3,0);
      m_empty=true;
    }; // default with type
    BoundarySpacePoint( BoundaryEnd_t type, const std::vector<float>& pos, const std::vector<float>& dir, const larcv::ImageMeta& );
    BoundarySpacePoint( BoundaryEnd_t type, const std::vector<float>& pos, const larcv::ImageMeta& );
    

#ifndef __CINT__
#ifndef __CLING__
    // move constructor
    BoundarySpacePoint( BoundaryEnd_t type, std::vector<BoundaryEndPt>&& endpts, const larcv::ImageMeta& meta  )  // type and endpt vector
      : std::vector<BoundaryEndPt>(std::move(endpts)) {
      boundary_type = type;
      m_empty = false;
      setup(meta);
    };
  
  // move constructor
  BoundarySpacePoint( BoundaryEnd_t type, std::vector<BoundaryEndPt>&& endpts, float x, float y, float z  )  // type and endpt vector and 3D poosition
    : std::vector<BoundaryEndPt>(std::move(endpts)) {
    boundary_type = type;
    m_pos.resize(3,0.0);
    m_dir.resize(3,0.0);
    m_pos[0] = x;
    m_pos[1] = y;
    m_pos[2] = z;
    m_empty = false;
  };
#endif
#endif

  // Copy Constructor
  //BoundarySpacePoint( const BoundarySpacePoint& src );
  
  virtual ~BoundarySpacePoint() {};

    BoundaryEnd_t type() const { return boundary_type; }
    const std::vector<float>& pos() const { return m_pos; }
    const std::vector<float>& dir() const { return m_dir; }
    void pos( std::vector<float>& pos_ ) {
      if ( pos_.size()!=3 ) {
    	throw std::runtime_error("BoundarySpacePoint::pos( vector )[error] Position vector must be 3-D (xyz in detector coordinates)");
      }
      m_pos = pos_;
      m_empty = false;
    }
    void dir( std::vector<float>& dir_ ) {
      if ( dir_.size()!=3 ) {
    	throw std::runtime_error("BoundarySpacePoint::dir( vector )[error] Position vector must be 3-D (xyz in detector coordinates)");
      }
      m_dir = dir_;
    }
    void setZY( float z, float y ) { m_pos[1] = y; m_pos[2] = z; m_empty = false; }
    void setX( float x ) { m_pos[0] = x; };
    ///< hmm, might be better to infer ZY position from endpt vector at constructor
    bool isempty() const { return m_empty; }
    float dwall() const;
    bool operator==( const BoundarySpacePoint& rhs ) const {
      if ( rhs.size()!=size() )
        return false;
      for (size_t p=0; p<(*this).size(); p++) {
        if ( (*this).at(p).col!=rhs.at(p).col || (*this).at(p).row!=rhs.at(p).row ) {
          return false;
        }
      }
      return true;
    };

    int tick( const larcv::ImageMeta& meta ) const;
    std::vector<int> wires( const larcv::ImageMeta& meta ) const;
    std::string printImageCoords( const larcv::ImageMeta& meta ) const;
    void setFlashIndex( int ivec, int idx, larlite::opflash* popfl=NULL );
    const BoundaryFlashIndex& getFlashIndex() const { return m_flashidx; };

protected:
    BoundaryEnd_t boundary_type;
    std::vector<float> m_pos;
    std::vector<float> m_dir;
    void setup( const larcv::ImageMeta& meta );
    bool m_empty; ///< flag that marks that the position has not been filled. used to indicate error states sometimes.
    BoundaryFlashIndex m_flashidx;
  };


}


#endif

