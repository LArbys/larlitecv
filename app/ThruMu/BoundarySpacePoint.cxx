#include "BoundarySpacePoint.h"
#include <cmath>

namespace larlitecv {

	float BoundarySpacePoint::dwall() const {
		float dwall = 0.0;
		float dy,dz;
		switch( type() ) {
			case larlitecv::kTop:
				dwall = 118.0-m_pos[1];
				break;
			case larlitecv::kBottom:
				dwall = 118.0+m_pos[1];
				break;
			case larlitecv::kUpstream:
				dwall = m_pos[2];
				break;
			case larlitecv::kDownstream:
				dwall = 1037-m_pos[2];
				break;
			case larlitecv::kAnode:
			case larlitecv::kCathode:
			case larlitecv::kImageEnd:
				dy = ( fabs(118.0-m_pos[1])<fabs(118.0+m_pos[1]) ) ? 118.0-m_pos[1] : 118.0+m_pos[1];
				dz = ( fabs(m_pos[2]) < fabs(1037-m_pos[2]) ) ? m_pos[2] : 1037-m_pos[2];
				dwall = ( fabs(dy)<fabs(dz ) ) ? dy : dz;
				break;
			default:
				std::runtime_error("BoundarySpacePoint::dwall[error] cannot calculate dwall for undefined boundary type");
				break;
		}
		return dwall;
	}

}