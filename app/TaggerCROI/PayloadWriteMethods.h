#ifndef __PAYLOAD_WRITE_METHODS_H__
#define __PAYLOAD_WRITE_METHODS_H__

/* --------------------------------------------------
// Collection of Functions for writing payload data 
// to disk via larlitecv::DataCoordinator
//
// ------------------------------------------------*/

// larlitecv
#include "Base/DataCoordinator.h"
#include "TaggerCROITypes.h"
#include "TaggerCROIAlgoConfig.h"

namespace larlitecv {

	void WriteInputPayload( const InputPayload&, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco );

	void WriteThruMuPayload( const ThruMuPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco );	

	void WriteStopMuPayload( const StopMuPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco );		

	void WriteCROIPayload( const CROIPayload& data, const InputPayload& inputdata, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco );	

}


#endif
