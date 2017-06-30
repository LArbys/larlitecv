#include "BoundaryMuonTaggerTypes.h"


namespace larlitecv {

  std::string BoundaryEndNames( const BoundaryEnd_t endpt )
  {
    switch (endpt) {
    case kTop:
      return "top";
    case kBottom:
      return "bottom";
    case kUpstream:
      return "upstream";
    case kDownstream:
      return "downstream";
    case kAnode:
      return "anode";
    case kCathode:
      return "cathode";
    case kImageEnd:
      return "image_end";
    default:
      return "unknown";
    }
    return "unknown";
  }
  
}
