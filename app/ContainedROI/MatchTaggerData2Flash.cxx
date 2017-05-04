#include "MatchTaggerData2Flash.h"

// larlite
#include "LArUtil/LArProperties.h"

namespace larlitecv {

  void MatchTaggerData2Flash( std::vector< TaggerFlashMatchData >& taggerdata_v, const std::vector< larlite::event_opflash* >& opflash_v,
			      const std::vector< BoundarySpacePoint >& anode_spacepts, const std::vector< BoundarySpacePoint >& cathode_spacepts,
			      const float max_dist ) {
    // In the case that we make 3D clusters and lose track of the boundary crossing types, we use this function to re-associate it for anode/cathode crossings
    // This happens when we recluster the tracks
    // This function also associates the opflashes with the track clusters as well.
    // We'll use this for CROI selection.

    enum MatchType_t { kAnode=0, kCathode, kNope };
    std::vector< const std::vector< BoundarySpacePoint >* > spacept_list(2,0);
    spacept_list[0] = &anode_spacepts;
    spacept_list[1] = &cathode_spacepts;

    const float usec_per_tick = 0.5;
    const float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*usec_per_tick;
    float dtick_drift = 258.0/0.111436/0.5 - 75.0;

    std::vector< int > startpt_type( taggerdata_v.size(), kNope );
    std::vector< int > endpt_type( taggerdata_v.size(), kNope );    

    // MATCH END POINTS TO BOUNDARY SPACE POINTS
    for ( size_t itagger=0; itagger<taggerdata_v.size(); itagger++ ) {
      auto& taggercluster = taggerdata_v[itagger];
      
      // first match end points to the space pts
      int matchtype = kAnode;
      float bestdist[2] = { 1.0e6, 1.0e6 };
      for ( auto const& pspacept_v : spacept_list ) {
	for ( auto const& spt : *pspacept_v ) {
	  
	  float dist[2] = {0.}; // [start,end]
	  for (int i=0; i<3; i++) {
	    dist[0] += ( spt.pos()[i] - taggercluster.m_track3d.Vertex()[i] )*( spt.pos()[i] - taggercluster.m_track3d.Vertex()[i] );
	    dist[1] += ( spt.pos()[i] - taggercluster.m_track3d.End()[i] )*( spt.pos()[i] - taggercluster.m_track3d.End()[i] );
	  }
	  dist[0] = sqrt(dist[0]);
	  dist[1] = sqrt(dist[1]);
	  
	  if ( dist[0]<max_dist ) {
	    if ( dist[0]<bestdist[0] ) {
	      startpt_type[itagger] = matchtype;
	      bestdist[0] = dist[0];
	    }
	  }
	  if ( dist[1]<max_dist ) {
	    if ( dist[1]<bestdist[1] ) {
	      endpt_type[itagger] = matchtype;
	      bestdist[1] = dist[1];
	    }
	  }
	  
	}//end of spacept vector loop
	matchtype++;
      }//end of spacept list loop
    }//end of tagger data loop

	// std::cout << "  [Start Flash]:"
	// 	  << " tick=" << 3200.0+taggerdata.m_track3d.Vertex()[0]/cm_per_tick
	// 	  << " start-z=" << taggerdata.m_track3d.Vertex()[2];
	// if ( !taggerdata.hasStartFlash() )
	//   std::cout << " no flash ";
	// else {
	//   std::cout << " tick=" << 3200.0+taggerdata.m_pstart_flash->Time()/0.5 + 15 << "/"
	// 	    << 3200.0+taggerdata.m_pstart_flash->Time()/0.5 + dtick_drift
	// 	    << " pe=" << taggerdata.m_pstart_flash->TotalPE()
	// 	    << " z=" << taggerdata.m_pstart_flash->ZCenter();
	// }
	
	// std::cout << "  [End Flash]: "
	// 	  << " tick=" << 3200.0+taggerdata.m_track3d.End()[0]/cm_per_tick
	// 	  << " end-z=" << taggerdata.m_track3d.End()[2];
	// if ( !taggerdata.hasEndFlash() )
	//   std::cout << " none ";
	// else {
	//   std::cout << " tick=" << 3200.0+taggerdata.m_pend_flash->Time()/0.5 + 15 << "/"
	// 	    << 3200.0+taggerdata.m_pend_flash->Time()/0.5 + dtick_drift
	// 	    << " pe=" << taggerdata.m_pend_flash->TotalPE()
	// 	    << " z=" << taggerdata.m_pend_flash->ZCenter();
	// }
	// std::cout << std::endl;


    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << "Rematched end point and flash data to reclustered tracks" << std::endl;
    
    // NOW ASSOCIATE THE FLASH
    for ( size_t itagger=0; itagger<taggerdata_v.size(); itagger++ ) {
      auto& taggercluster = taggerdata_v[itagger];
      std::cout << "  [" << itagger << "] ";
      //  START POINT
      if ( startpt_type[itagger]!=kNope ) {
	// try to find the flash that made this tag
	float tick = taggercluster.m_track3d.Vertex()[0]/cm_per_tick + 3200.0;
	const larlite::opflash* pbest_flash = NULL;
	float dtick = -1;
	float flash_tick = -1;
	if ( startpt_type[itagger]==kAnode ) {
	  std::cout << "[ Start (" << taggercluster.m_track3d.Vertex()[0] << "," << taggercluster.m_track3d.Vertex()[2] << "): Anode tick=" << tick;
	  // This start point matches an anode space point
	  for ( auto const& pflash_v : opflash_v ) {
	    for ( auto const& flash : *pflash_v ) {
	      float flash_tick_target = 3200.0 + flash.Time()/usec_per_tick + 15; // has to match ThruMu/FlashMuonTaggerAlgo (fix this hard coding)
	      if ( pbest_flash==NULL || dtick>fabs( flash_tick_target-tick ) ) {
		dtick = fabs( flash_tick_target-tick );
		pbest_flash = &flash;
		flash_tick = flash_tick_target;
	      }
	    }
	  }
	}
	else if ( startpt_type[itagger]==kCathode ) {
	  std::cout << "[ Start (" << taggercluster.m_track3d.Vertex()[0] << "," << taggercluster.m_track3d.Vertex()[2] << "): Cathode tick=" << tick;	  
	  // This start point matches an anode space point
	  for ( auto const& pflash_v : opflash_v ) {
	    for ( auto const& flash : *pflash_v ) {
	      float flash_tick_target = 3200.0 + flash.Time()/usec_per_tick + dtick_drift; // has to match ThruMu/FlashMuonTaggerAlgo (fix this hard coding)
	      if ( pbest_flash==NULL || dtick>fabs( flash_tick_target-tick ) ) {
		dtick = fabs( flash_tick_target-tick );
		pbest_flash = &flash;
		flash_tick = flash_tick_target;
	      }
	    }
	  }
	}
	
	if ( pbest_flash ) {
	  taggercluster.setStartFlash( pbest_flash );
	  std::cout << " flashtick=" << flash_tick << " dtick=" << dtick << " z=" << pbest_flash->ZCenter() << " pe=" << pbest_flash->TotalPE();
	}
      }// end of start point if
      else {
	std::cout << "[ Start: NoFlash";
      }

      std::cout << " ] ";
      
      // END POINT
      if ( endpt_type[itagger]!=kNope ) {
	// try to find the flash that made this tag
	float tick = taggercluster.m_track3d.End()[0]/cm_per_tick + 3200.0;
	const larlite::opflash* pbest_flash = NULL;	
	float dtick = -1;
	float flash_tick = -1;
	if ( endpt_type[itagger]==kAnode ) {
	  std::cout << "[ End(" << taggercluster.m_track3d.End()[0] << "," << taggercluster.m_track3d.End()[2] << "): Anode tick=" << tick;
	  // This start point matches an anode space point
	  for ( auto const& pflash_v : opflash_v ) {
	    for ( auto const& flash : *pflash_v ) {
	      float flash_tick_target = 3200.0 + flash.Time()/usec_per_tick + 15; // has to match ThruMu/FlashMuonTaggerAlgo (fix this hard coding)
	      if ( pbest_flash==NULL || dtick>fabs( flash_tick_target-tick ) ) {
		dtick = fabs( flash_tick_target-tick );
		pbest_flash = &flash;
		flash_tick = flash_tick_target;
	      }
	    }
	  }
	}
	else if ( endpt_type[itagger]==kCathode ) {
	  std::cout << "[ End(" << taggercluster.m_track3d.End()[0] << "," << taggercluster.m_track3d.End()[2] << "): Cathode tick=" << tick;	  
	  // This start point matches an anode space point
	  for ( auto const& pflash_v : opflash_v ) {
	    for ( auto const& flash : *pflash_v ) {
	      float flash_tick_target = 3200.0 + flash.Time()/usec_per_tick + dtick_drift; // has to match ThruMu/FlashMuonTaggerAlgo (fix this hard coding)
	      if ( pbest_flash==NULL || dtick>fabs( flash_tick_target-tick ) ) {
		dtick = fabs( flash_tick_target-tick );
		pbest_flash = &flash;
		flash_tick = flash_tick_target;		
	      }
	    }
	  }
	}
	if ( pbest_flash ) {
	  taggercluster.setEndFlash( pbest_flash );
	  std::cout << " flashtick=" << flash_tick << " dtick=" << dtick << " z=" << pbest_flash->ZCenter() << " pe=" << pbest_flash->TotalPE();	  
	}	
      }//end of endpt block
      else {
	std::cout << "[ End: NoFlash";
      }

      std::cout << " ]" << std::endl;
      
    }//end of loop over tagger data
    std::cout << "----------------------------------------------------------" << std::endl;    
    
  }//end of function

}
