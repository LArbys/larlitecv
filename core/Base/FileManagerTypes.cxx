#include "FileManagerTypes.h"

namespace larlitecv {

  std::ostream& operator<<(std::ostream &os, RSE const& m ) {
    std::stringstream ss;
    ss << "RSE(" << m.run << ", " << m.subrun << ", " << m.event;
    if ( m.subevent!=0 )
    	ss << ", " << m.subevent << ")";
    else
    	ss << ")";
    return (os << ss.str());
  };

}
