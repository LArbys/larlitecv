#include "PayloadWriteMethods.h"

#include <vector>

// larlite
#include "DataFormat/opflash.h"

// larcv
#include "DataFormat/DataFormatTypes.h"
#include "DataFormat/EventImage2D.h"

namespace larlitecv {

	void WriteInputPayload( const InputPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

		// input
		if ( config.input_write_cfg.get<bool>("SaveInputTPC") ) {
  		larcv::EventImage2D* evout_img = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "modimg" );
	  	for ( auto const& img : data.img_v )
		  	evout_img->Append( img );
		}

		// opflash
		if ( config.input_write_cfg.get<bool>("SaveOpflash") ) {
  		larlite::event_opflash* evout_opflash = (larlite::event_opflash*)dataco.get_larlite_data( larlite::data::kOpFlash, "opflash" );
  		for ( auto const& opflash_v : data.opflashes_v ) {
  			for ( auto const& opflash : *opflash_v )
    			evout_opflash->push_back( opflash );
  		}
		}

		// badch
		if ( config.input_write_cfg.get<bool>("SaveBadChImage") ) {
			larcv::EventImage2D* evout_badch = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "badch" );
			for ( auto const& img : data.badch_v )
				evout_badch->Append( img );
		}

		// gapchs
		if ( config.input_write_cfg.get<bool>("SaveGapChImage") ) {
			larcv::EventImage2D* evout_gapchs = (larcv::EventImage2D*)dataco.get_larcv_data( larcv::kProductImage2D, "gapchs" );
			for ( auto const& img : data.gapch_v )
				evout_gapchs->Append( img );
		}

	}


	void WriteThruMuPayload( const ThruMuPayload& data, const TaggerCROIAlgoConfig& config, DataCoordinator& dataco ) {

	}
}