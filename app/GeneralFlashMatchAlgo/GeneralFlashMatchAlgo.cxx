// (All of the functions for flash matching from 'ContainedROI/TaggerFlashMatchAlgo.h' are included in "GeneralFlashMatchAlgo.h".
#include "GeneralFlashMatchAlgo.h"

namespace larlitecv {

  GeneralFlashMatchAlgo::GeneralFlashMatchAlgo()
    : m_pmtweights("geoinfo.root")
  {}

  GeneralFlashMatchAlgo::GeneralFlashMatchAlgo( GeneralFlashMatchAlgoConfig& config )
    : m_config(config), m_pmtweights("geoinfo.root") {
    
    setVerbosity( m_config.verbosity );
    m_flash_matcher.Configure( m_config.m_flashmatch_config );
    
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();
    for ( size_t opch=0; opch<32; opch++) {
      m_opch_from_opdet.insert( std::make_pair<int,int>( geo->OpDetFromOpChannel(opch), opch ) );
    }
  }
  
  // A function that will find the flash indices that correspond to flashes that do not correspond to anode/cathode piercing endpoints.
  // Inputs: opflash_v: all of the opflash objects from the event.
  //         anode_flash_idx_v: the indices of the flashes that correspond to boundary points at the anode.
  //         cathode_flash_idx_v: the indices of the flashes that correspond to boundary points at the cathode.
  std::vector < int > GeneralFlashMatchAlgo::non_anodecathode_flash_idx( const std::vector < larlite::opflash* >& opflash_v,
									 std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v ) {
    // Declare a vector for the flashes that are neither anode or cathode crossing.
    std::vector < int > non_ac_flash_idx;
    
    // Loop through the opflash information from the event in order to separate the anode/cathode crossing flashes from the others.
    for (size_t non_ac_iter = 0; non_ac_iter < opflash_v.size(); non_ac_iter++) {
      
      bool belongs_to_ac_flash = false;
      
      // Loop through the anode flashes.
      for (size_t anode_flsh_i = 0; anode_flsh_i < anode_flash_idx_v.size(); anode_flsh_i++) {

	// Compare the two flash indices.
	if ((int)non_ac_iter == anode_flash_idx_v[anode_flsh_i]) belongs_to_ac_flash = true;

      }

      // Loop through the cathode flashes.
      for (size_t cathode_flsh_i = 0; cathode_flsh_i < cathode_flash_idx_v.size(); cathode_flsh_i++) {

        // Compare the two flash indices.
	if ((int)non_ac_iter == cathode_flash_idx_v[cathode_flsh_i]) belongs_to_ac_flash = true;

      }

      // If it was not the same index value as any of the anode/cathode flash indices, then you can append it onto the 'non_ac_flash_idx' vector.
      if (!belongs_to_ac_flash) non_ac_flash_idx.push_back(non_ac_iter);

    }

    // Return the vector of non anode-piercing/cathode-piercing track endpoints.
    return non_ac_flash_idx;

  }


  // A function that will generate a list of opflashes that correspond to the list of indices that are put in.
  // Inputs: full_opflash_v: This is the full list of opflashes from which we are trying to filter
  //   opflashes of the same denomination,
  // i.e. all anode-piercing flashes, all cathode-piercing flashes, all non-side-piercing flashes, etc.
  //         idx_v: This is the list of indices of the flashes that you are interested in.
  // You want to place the flashes located in the 'full_opflash_v' vector at each index in this list
  // in the output list, which contains the opflashes that we are interested in.
  std::vector <larlite::opflash*> GeneralFlashMatchAlgo::generate_single_denomination_flash_list( std::vector< larlite::opflash > full_opflash_v, std::vector < int > idx_v) {

    // Declare a list of flashes that will be the output, which is the list of the certain denomination of flashes.
    std::vector < larlite::opflash* > single_denomination_flash_list;
    single_denomination_flash_list.clear();
    
    for (size_t iflash = 0; iflash < idx_v.size(); iflash++) {
      larlite::opflash* opflash = &(full_opflash_v.at(idx_v.at(iflash)));
      single_denomination_flash_list.push_back( opflash );
    }

    // Return this list of the opflash information of one denomination.
    return single_denomination_flash_list;

  }


  // Define a function that will flashmatch the contained ROIs in the TaggerCROI selection.
  // Inputs: 'flashdata_v' - This is of type 'larlitecv::TaggerFlashMatchData' and contains the cluster type (ThruMu, StopMu, Contained/Untagged), the 'Pixel2DCluster', and the larlite::track version of each object in the TPC.
  //         'opflashes_v' - This is a vector of type 'event_opflash', which itself is a vector of type 'opflash'.
  //         'flashdata_selected_v' - This is a vector of type 'int'.

  std::vector<larcv::ROI> GeneralFlashMatchAlgo::FindFlashMatchedContainedROIs( const std::vector< larlitecv::TaggerFlashMatchData>& flashdata_v, const std::vector< larlite::event_opflash* >& opflashes_v, std::vector< int >& flashdata_selected_v ) {

    // the output                                                                                                                                                                                      
    std::vector<larcv::ROI> roi_v;

    // clear state                                                                                                                                                                                     
    if (flashdata_selected_v.size()!=flashdata_v.size()) {
      flashdata_selected_v.resize( flashdata_v.size(), 0 );
    }

    // get the in-beam flashes                                                                                                                                                                         
    std::vector<flashana::Flash_t> data_flashana = GetInTimeFlashana( opflashes_v );
    if ( data_flashana.size()==0 ) {
      return roi_v;
    }
   
    // choose contained candidates
    m_passes_containment.resize( flashdata_v.size(), 0 );
    m_passes_flashmatch.resize( flashdata_v.size(), 0 );
    m_passes_totpe.resize( flashdata_v.size(), 0 );
    m_passes_cosmicflash_ratio.resize( flashdata_v.size(), 0 );

    m_min_chi2.clear();
    m_opflash_hypos.clear();
    m_totpe_peratio.clear();
    m_cosmicflash_ratio_dchi.clear();
    m_min_chi2.reserve( flashdata_v.size() );
    m_totpe_peratio.reserve( flashdata_v.size() );
    m_cosmicflash_ratio_dchi.reserve( flashdata_v.size() );

    for ( size_t i=0; i<flashdata_v.size(); i++ ) {

      if ( m_verbosity>=2 ) {
	std::cout << " Candidate #" << i << ", ";
	if ( flashdata_v.at(i).m_type==TaggerFlashMatchData::kThruMu )
	  std::cout << "ThruMu";
	else if ( flashdata_v.at(i).m_type==TaggerFlashMatchData::kStopMu )
	  std::cout << "StopMu";
	else if ( flashdata_v.at(i).m_type==TaggerFlashMatchData::kUntagged )
	  std::cout << "Untagged";
	std::cout << ": ";
      }

      m_passes_containment[i] = ( IsClusterContained( flashdata_v.at(i) ) ) ? 1 : 0;
      
      if ( m_verbosity>=2 ) {
	if ( m_passes_containment[i] )
	  std::cout << " {contained}" << std::endl;
	else
	  std::cout << "{uncontained}" << std::endl;
      }

      // flash-match candidates 
      // Declare a qcluster here that will be filled in the 'ExpandQClusterStartingWithLarliteTrack' function.
      flashana::QCluster_t qcluster;

      // We do not want to extend this qcluster outside the TPC, both 'extend_start' and 'extend_end' are set to false.
      // Ask the 'flashdata_v' object for 'm_track3d', which is the 3D larlite track of the 'TaggerFlashMatchData' object.
      ExpandQClusterStartingWithLarliteTrack(qcluster, flashdata_v.at( i ).m_track3d, 0.0, false, false);

      // in-time flash matching                                                                                                                                                                      
      float totpe_data = 0;
      float totpe_hypo = 0;
      m_passes_flashmatch[i] = ( DoesQClusterMatchInTimeFlash( data_flashana, qcluster, totpe_data, totpe_hypo )  ) ? 1 : 0;

      // totpe                                                                                                                                                                                        
      m_passes_totpe[i] = ( DoesTotalPEMatch( totpe_data, totpe_hypo ) ) ? 1 : 0;

      // cosmic versus in-time flash log-likelihood ratio                                                                                                                                             
      float dchi2=0;
      m_passes_cosmicflash_ratio[i] = ( DoesQClusterMatchInTimeBetterThanCosmic( data_flashana, qcluster, flashdata_v.at( i ), i, dchi2 ) ) ? 1 : 0;

      if ( m_verbosity>=2 ) {
	std::cout << "  ";

        if ( m_passes_containment[i] )
	  std::cout << " {contained}";
        else
	  std::cout << "{uncontained}";

        if ( m_passes_flashmatch[i]==1 ) {
	  std::cout << " {in-time flash-matched}";
        }
        else {
	  std::cout << " {not-matched}";
        }

        if ( m_passes_totpe[i]==1 )
	  std::cout << " {totpe matched}";
        else
	  std::cout << " {totpe fails}";

        if ( !m_passes_cosmicflash_ratio[i] )
	  std::cout << " {fails cosmic-flash match: dchi2=" << dchi2 << "}";
        else if ( dchi2!=0 )
	  std::cout << " {passes cosmic-flash match: dchi2=" << dchi2 << "}";
        else
	  std::cout << " { no flash end }";

	if ( m_passes_flashmatch[i] && m_passes_containment[i] && m_passes_cosmicflash_ratio[i] && m_passes_totpe[i] )
	  std::cout << " **PASSES**";
	std::cout << std::endl;
      }
    }//end of input data loop                                                                                                                                                                               

    for ( size_t i=0; i<flashdata_v.size(); i++) {
      if ( m_passes_containment[i] && m_passes_flashmatch[i] && m_passes_cosmicflash_ratio[i] && m_passes_totpe[i] ) {
        // Make ROI                                                                                                                                                                                         
        flashdata_selected_v[i] = 1;
      }
    }

    return roi_v;
  }



  std::vector<larlite::opflash> GeneralFlashMatchAlgo::SelectInTimeOpFlashes( const std::vector<larlite::event_opflash*>& opflashes_v ) {

    std::vector<larlite::opflash> beam_flashes;
    for ( auto const& ptr_ev_flash : opflashes_v ) {
      for ( auto const& opflash : *ptr_ev_flash ) {
        if ( opflash.TotalPE()<m_config.flashpe_thresh )
          continue;
        int tick = opflash.Time()/m_config.us_per_tick;
        if ( tick>=m_config.beam_tick_range[0] && tick <=m_config.beam_tick_range[1] ) {
          if ( m_verbosity>0 )
	    std::cout << "In-time flash found: " << opflash.Time() << "us from trigger. Tick=" << tick << std::endl;
          beam_flashes.push_back( opflash ); // our own copy                                                                                                                                           
	}
        else {
          if ( m_verbosity>0)
	    std::cout << "Rejected flash: " << opflash.Time() << "us from trigger. Tick=" << tick << std::endl;
        }
      }
    }

    if ( m_verbosity>0 )
      std::cout << "TaggerFlashMatchAlgo::SelectInTimeFlashes. Found " << beam_flashes.size() << " flashes." << std::endl;

    return beam_flashes;
  }



  bool GeneralFlashMatchAlgo::DoesTotalPEMatch( float totpe_data, float totpe_hypo )  {
    float frac_diff = fabs(totpe_data-totpe_hypo)/totpe_hypo;
    m_totpe_peratio.push_back(frac_diff);
    if ( frac_diff > m_config.totpe_sigma_cut )
      return false;
    return true;
  }

  bool GeneralFlashMatchAlgo::IsClusterContained( const TaggerFlashMatchData& data ) {
    // simple selection                                                                                                                                                                               
    // we make a cut on fiducial volume as a function of x                                                                                                                                             
    
    // first we need the bounding box of the points                                                                                                                                                    
    float bb[3][2] = {0}; // extremes in each dimension                                                                                                                                                
    float extrema[3][2][3] = {-1.0e6 }; // each extrema point stored here. (dim,min/max,xyz)                                                                                                           
    for ( int v=0; v<3; v++) {
      bb[v][0] = 1.0e6;
      bb[v][1] = -1.0e6;
    }

    const larlite::track& track = data.m_track3d;
    for ( size_t i=0; i<track.NumberTrajectoryPoints(); i++ ) {
      const TVector3& xyz = track.LocationAtPoint(i);
      for (int v=0; v<3; v++) {
        // minvalue                                                                                                                                                                                    
	if ( bb[v][0]>xyz[v] ) {
          bb[v][0] = xyz[v];
          for (int j=0; j<3; j++)
            extrema[v][0][j] = xyz[j];
        }
        // maxvalue                                                                                                                                                                                      
	if ( bb[v][1]<xyz[v] ) {
          bb[v][1] = xyz[v];
          for (int j=0; j<3; j++)
            extrema[v][1][j] = xyz[j];
        }
      }
    }

    // we adjust for space charge effects                                                                                                                                                              
    //std::vector<double> delta_ymin = m_sce.GetPosOffsets( extrema[1][0][0], -118.0, extrema[1][0][2] );                                                                                               
    //std::vector<double> delta_ymax = m_sce.GetPosOffsets( extrema[1][1][0],  118.0, extrema[1][1][2] );                                                                                               
    //std::vector<double> delta_zmin = m_sce.GetPosOffsets( extrema[2][0][0], extrema[2][0][1], 0.0    );                                                                                              
    //std::vector<double> delta_zmax = m_sce.GetPosOffsets( extrema[2][1][0], extrema[2][1][1], 1037.0 );                                                                                               
    std::vector<double> delta_ymin(3,0);
    std::vector<double> delta_ymax(3,0);
    std::vector<double> delta_zmin(3,0);
    std::vector<double> delta_zmax(3,0);


    if ( m_verbosity>=2 ) {
      std::cout << "bounds: x=[" << bb[0][0] << "," << bb[0][1] << "] "
                << " y=[" << bb[1][0] << "," << bb[1][1] << "] "
                << " z=[" << bb[2][0] << "," << bb[2][1] << "] "
                << " dy=[" << delta_ymin[1] << "," << delta_ymax[1] << "] "
                << " dz=[" << delta_zmin[2] << "," << delta_zmax[2] << "] ";
    }

    // x extrema, not a function of the other dimensions                                                                                                                                               
    if ( bb[0][0]<m_config.FVCutX[0] || bb[0][1]>m_config.FVCutX[1] )
      return false;
    if ( (bb[1][0]-delta_ymin[1]) < m_config.FVCutY[0] || (bb[1][1]-delta_ymax[1]) > m_config.FVCutY[1] )
      return false;
    if ( (bb[2][0]-delta_zmin[2]) < m_config.FVCutZ[0] || (bb[2][1]-delta_zmax[2]) > m_config.FVCutZ[1] )
      return false;

    // we made it!                                                                                                                                                                                     
    return true;
  }




  // Define a function that will return the in-time flashes from the event.
  std::vector<flashana::Flash_t> GeneralFlashMatchAlgo::GetInTimeFlashana( const std::vector<larlite::event_opflash*>& opflashes_v ) {

    std::vector<larlite::opflash> intime_opflashes_v = SelectInTimeOpFlashes( opflashes_v );

    std::vector<flashana::Flash_t> intime_flashana_v = MakeDataFlashes( intime_opflashes_v );

    return intime_flashana_v;
  }



  // Defining a function for if the qcluster matches the in-time flash.
  bool GeneralFlashMatchAlgo::DoesQClusterMatchInTimeFlash( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster, float& totpe_data, float& totpe_hypo ) {

    // Perform a loop to find the minimum chi2 match between the qcluster and any of the intime flashes.                                                                                              
    // Minimize this variable in the loop.                                                                                                                                                             
    // Declare the final values that these variables will hold.                                                                                                                                        
    float chi2                 = -1.0;
    totpe_data                 = 0.0;
    totpe_hypo                 = 0.0;

    // Declare the initial values that these variables will hold while it is being minimized.                                                                                                           
    float min_chi2_candidate   = 0.0;
    float totpe_data_candidate = 0.0;
    float totpe_hypo_candidate = 0.0;

    for ( auto& intime_flash: intime_flashes_v ) {

      // Convert the intime flash to an opflash object.
      // It will be converted back into a 'flashana::Flash_t' object in the other function, but this is just for compatability.
      larlite::opflash intime_opflash = MakeOpFlashFromFlash( intime_flash );

      min_chi2_candidate = generate_chi2_in_track_flash_comparison( qcluster, intime_opflash, totpe_data_candidate, totpe_hypo_candidate, 0 );

      if ( min_chi2_candidate < chi2 || chi2 < 0.0) {
        chi2       = min_chi2_candidate;
        totpe_data = totpe_data_candidate;
        totpe_hypo = totpe_hypo_candidate;
      }

    }


    // our cut needs to be position dependent, so we need a position. we get the q-weighted mean and the intervals.                                                                                   
    float totw = 0.;
    float mean[3] = {0};
    float min_x = 1.0e6;
    float max_x = -1.0e6;
    for ( auto const& qpt : qcluster ) {
      mean[0] += qpt.x*qpt.q;
      mean[1] += qpt.y*qpt.q;
      mean[2] += qpt.z*qpt.q;
      totw += qpt.q;
      if ( qpt.x < min_x ) min_x = qpt.x;
      if ( qpt.x > max_x ) max_x = qpt.x;
    }
    if ( totw==0 ) {
      return false;
    }
    else {
      for (int i=0; i<3; i++) mean[i] /= totw;
    }

    //  still a dumb cut                                                                                                                                                                               
    if ( mean[0] < 100.0 && chi2< m_config.flashmatch_chi2_cut*1.5 )
      return true;
    else if ( mean[0]>100.0 && chi2<m_config.flashmatch_chi2_cut )
      return true;
    else
      return false;

  }

  
  // Add a function that checks if a qcluster matches an in-time flash better than a cosmic.

  bool GeneralFlashMatchAlgo::DoesQClusterMatchInTimeBetterThanCosmic( const std::vector<flashana::Flash_t>& intime_flashes_v, const flashana::QCluster_t& qcluster,
                                                                      const larlitecv::TaggerFlashMatchData& taggertrack, const int trackidx, float& dchi2 ) {
    if ( !taggertrack.hasStartFlash() && !taggertrack.hasEndFlash() ) {
      // No flash associated to this track.                                                                                                                                                            
      return true;
    }

    // check that the matching flash is not the intime flash                                                                                                                                           
    float intime_min_chi2 = m_min_chi2[trackidx]; // we cheat and use this stored value instead of recalculating                                                                                      

    float totpe_cosmicflash = 0.;
    float totpe_hypoflash = 0;

    float cosmicflash_start_chi2 = 1.0e6;
    if ( taggertrack.hasStartFlash() ) {
      std::vector< larlite::opflash > cosmicflash_start_v;
      float start_tick = taggertrack.m_pstart_flash->Time()/m_config.us_per_tick;
      float start_ly = 0;
      if ( start_tick<m_config.beam_tick_range[0] || start_tick>m_config.beam_tick_range[1] )
        start_ly = 1; // Use cosmic disc.  This is used in the context now as a constant multiplied by the same light yield.
      cosmicflash_start_v.push_back( *(taggertrack.m_pstart_flash ) );
      cosmicflash_start_chi2 = generate_chi2_in_track_flash_comparison( qcluster, cosmicflash_start_v.at(0), totpe_cosmicflash, totpe_hypoflash, start_ly ); 
    }

    float cosmicflash_end_chi2 = 1.0e6;
    if ( taggertrack.hasEndFlash() ) {
      std::vector< larlite::opflash > cosmicflash_end_v;
      float end_ly = 0;
      float end_tick = taggertrack.m_pend_flash->Time()/m_config.us_per_tick;
      if ( end_tick<m_config.beam_tick_range[0] || end_tick>m_config.beam_tick_range[1] )
        end_ly = 1; // use cosmic disc.  This is used in the context now as a constant multiplied by the same light yield. 
      cosmicflash_end_v.push_back( *(taggertrack.m_pend_flash ) );
      cosmicflash_end_chi2 = generate_chi2_in_track_flash_comparison( qcluster, cosmicflash_end_v.at(0), totpe_cosmicflash, totpe_hypoflash, end_ly );
    }

    float cosmicflash_chi2 = ( cosmicflash_start_chi2<cosmicflash_end_chi2 ) ? cosmicflash_start_chi2 : cosmicflash_end_chi2;

    dchi2 = intime_min_chi2 - cosmicflash_chi2;
    m_cosmicflash_ratio_dchi.push_back( dchi2 );

    if ( cosmicflash_chi2 < intime_min_chi2 ) {
      // matches cosmic flash better                                                                                                                                                                
      return false;
    }
    else {
      return true;
    }

    return true; // never gets here                                                                                                                                                                   
  }


  // Declare a function that will store all of the flashes for an event in a running list rather than separating them according to their flash producer.
  // This list assumes that there are two groups of flashes, one produced using 'simpleFlashBeam' and the second produced using 'simpleFlashCosmic', in that order.
  // Input list: 'opflashes_v', which contains the sets of opflash products in different vectors within the entire vector, with each corresponding to a different flash producer.
  std::vector< larlite::opflash > GeneralFlashMatchAlgo::generate_single_opflash_vector_for_event( std::vector< larlite::event_opflash* > opflashes_v ) {

    // Declare an output list of opflash products that will be filled with the information in the input vector, 'opflashes_v'.
    std::vector< larlite::opflash > single_opflash_vector;
    single_opflash_vector.clear();

    // Use an iterator to resize the output vector.
    int total_opflash_list_iter = 1;

    // Loop through the input vector and each vector within to generate 'single_opflash_vector' from its contents.
    for ( auto& ptr_event_flash : opflashes_v ) {

      // Generate an object for the opflash producer that we are currently looking at.
      auto& current_event_opflash_object = *ptr_event_flash;

      // Loop through the inner vector of each 'event_opflash*> producer contained within 'opflashes_v' to fill the output 'single_opflash_vector'.

      for ( size_t inner_opflash_iter = 0; inner_opflash_iter < current_event_opflash_object.size(); ++inner_opflash_iter ) {

	single_opflash_vector.resize( total_opflash_list_iter );

	// Append the opflash pointer located at this iteration in the inner loop over the 'opflashes_v' vector to the 'single_opflash_vector' output.
	larlite::opflash single_opflash                    = current_event_opflash_object.at( inner_opflash_iter );
	single_opflash_vector[total_opflash_list_iter - 1] = single_opflash;

	++total_opflash_list_iter;

      }

    }

    // Return this single denomination flash list.
    return single_opflash_vector;

  }

  // Declare a function that will return a list with our numbering scheme for the producer that made the flashes in the event.
  // This list assumes that there are two groups of flashes, one produced using 'simpleFlashBeam' and the second produced using 'simpleFlashCosmic', in that order.
  // This function will return a '0' for 'simpleFlashBeam' and a '1' for 'simpleFlashCosmic'.
  // Input list: 'opflashes_v', which contains the sets of opflash products in different vectors within the entire vector, with each corresponding to a different flash producer.
  std::vector< int > GeneralFlashMatchAlgo::generate_single_opflash_idx_vector_for_event( std::vector< larlite::event_opflash* > opflashes_v ) {

    // Declare an output list of opflash products that will be filled with the information for opflash producer, which takes the form of a vector of integer values.
    std::vector< int > single_opflash_idx_vector;
    single_opflash_idx_vector.clear();

    // Use an iterator to fill 'single_opflash_idx_vector' with.
    int single_opflash_idx = 0;

    // Start the outer loop.
    for ( auto& opflash_set: opflashes_v ) {

      // Reformat the function like this so that it compiles.
      auto& current_opflash_object = *opflash_set;

      // Start the inner loop.
      for ( auto& opflash: current_opflash_object ) {

	single_opflash_idx_vector.push_back( single_opflash_idx );

      }

      // Increment 'single_opflash_idx'.
      ++single_opflash_idx;

    }

    // Return 'single_opflash_idx_vector'.
    return single_opflash_idx_vector;

  }

  // A function that will generate tracks from the 'BMTrackCluster3D' objects that are formed from each pass of the tagger.
  // Inputs: trackcluster3d_v: A vector of 'BMTrackCluster3D' objects that were formed from a pass of the tagger.
  std::vector < larlite::track >  GeneralFlashMatchAlgo::generate_tracks_between_passes(const std::vector< larlitecv::BMTrackCluster3D > trackcluster3d_v) {

    // We will use the 'makeTrack()' functionality of these 'BMTrackCluster3D' objects to turn them into a vector of 'larlite::track' objects.
    std::vector < larlite::track > output_tracks;

    // Loop through the list of 'BMTrackCluster3D' objects and turn them into a vector 'larlite::track' objects.
    for (auto const& trackcluster3d : trackcluster3d_v) {

      output_tracks.push_back(trackcluster3d.makeTrack());

    }

    // Return this list of tracks.
    return output_tracks;

  }

  // Final function to extend the qcluster, analogous to the function of 'ExpandQCluster' found in 'MCQCluster'.  This function will take the entire larlite track.
  // Inputs: 'qcluster'      - an empty qcluster that is formed from a larlite track.
  //         'larlite_track' - the input track from which the qcluster is formed.
  //         'extension'     - the amount by which the cluster is extended outside the TPC (in cm).  This number is arbitrarily set to 10 m so that the track is extended to the edge of the cryostat and eyond.            
  //         'extend_start'  - a boolean to allow the user to refrain from extending the start of track if that is desired.
  //         'extend_end'    - a boolean to allow the user to refrain from extending the end of the track if that is desired.
  void GeneralFlashMatchAlgo::ExpandQClusterStartingWithLarliteTrack( flashana::QCluster_t& qcluster, const larlite::track& larlite_track, double extension, bool extend_start, bool extend_end ) {

    // Continue over the track if it has less than 2 trajectory points in its length.
    if ( larlite_track.NumberTrajectoryPoints() < 2 ) return;

    // Declare the two points that will be used in the function.
    ::geoalgo::Vector pt0(0.,0.,0.);
    ::geoalgo::Vector pt1(0.,0.,0.);

    // Declare an object for 'lightpath', which will extend the qcluster.
    flashana::LightPath lightpath;

    //
    // Add body.
    //
    for ( size_t step_idx = 0; step_idx < larlite_track.NumberTrajectoryPoints() - 1; ++step_idx ) {

      pt0[0] = larlite_track.LocationAtPoint(step_idx)[0];
      pt0[1] = larlite_track.LocationAtPoint(step_idx)[1];
      pt0[2] = larlite_track.LocationAtPoint(step_idx)[2];

      pt1[0] = larlite_track.LocationAtPoint(step_idx+1)[0];
      pt1[1] = larlite_track.LocationAtPoint(step_idx+1)[1];
      pt1[2] = larlite_track.LocationAtPoint(step_idx+1)[2];

      lightpath.QCluster( pt0, pt1, qcluster );

    }

    // Initialize the first point on the track's trajectory.
    pt0[0] = larlite_track.LocationAtPoint(0)[0]; pt0[1] = larlite_track.LocationAtPoint(0)[1]; pt0[2] = larlite_track.LocationAtPoint(0)[2];

    // Starting Point Analysis.
    // Check to see if this point is at the edge of the TPC volume.  I will here use a resolution value of 10.0 cm.
    if ( isNearActiveVolumeEdge(pt0, 10.0) && extend_start ) {
    
      int next_point_idx = 1;

      do {
	pt1[0] = larlite_track.LocationAtPoint( next_point_idx )[0]; pt1[1] = larlite_track.LocationAtPoint( next_point_idx )[1]; pt1[2] = larlite_track.LocationAtPoint( next_point_idx )[2];
	++next_point_idx;
      } while ( fabs( pt0[0] - pt1[0] ) < 0.001 && fabs( pt0[1] - pt1[1] ) < 0.001 && fabs( pt0[2] - pt1[2] ) < 0.001 && next_point_idx < int( larlite_track.NumberTrajectoryPoints() ) );

      // Continue on only if 'next_point_idx' is less than the number of trajectory points on the track.
      if ( next_point_idx < larlite_track.NumberTrajectoryPoints() ) {
	// Turn 'pt1' into a direction vector (not normalized).
	pt1 = pt0 - pt1;
	pt0    = pt0 + pt1.Dir()*extension;
	pt1[0] = larlite_track.LocationAtPoint(0)[0]; pt1[1] = larlite_track.LocationAtPoint(0)[1]; pt1[2] = larlite_track.LocationAtPoint(0)[2];

	// Extend the qcluster with 'lightpath' functionality.
	lightpath.QCluster( pt0, pt1, qcluster );

      }

  }
    
  // Initialize the last point on the track's trajectory.
  pt0[0] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[0]; pt0[1] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[1]; pt0[2] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[2];

  // Ending Point Analysis.                                                                                                                                                                               
  // Check to see if this point is on the edge of the TPC volume.  I will here use a resolution value of 10.0 cm.                                                                                         
  if ( isNearActiveVolumeEdge(pt0, 10.0) && extend_end ) {
    
    int next_point_idx = int(larlite_track.NumberTrajectoryPoints() - 2);
 
    do { 
      pt1[0] = larlite_track.LocationAtPoint( next_point_idx )[0]; pt1[1] = larlite_track.LocationAtPoint( next_point_idx )[1]; pt1[2] = larlite_track.LocationAtPoint( next_point_idx )[2];
    --next_point_idx;
    } while ( fabs( pt0[0] - pt1[0] ) < 0.001 && fabs( pt0[1] - pt1[1] ) < 0.001 && fabs( pt0[2] - pt1[2] ) < 0.001 && next_point_idx > -1 );

    // Only perform the next part of the loop if 'next_point_idx' is greater than -1.
    if ( next_point_idx > -1 ) {
      // Turn 'pt1' into a direction vector (not normalized).
      pt1 = pt0 - pt1;
      pt0 = pt0 + pt1.Dir()*extension;
      pt1[0] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[0]; pt1[1] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[1]; pt1[2] = larlite_track.LocationAtPoint( larlite_track.NumberTrajectoryPoints() - 1 )[2];

      // Extend the qcluster with the 'lightpath' functionality.                                                                                                                                      
      lightpath.QCluster( pt0, pt1, qcluster );
    }

  }

  return;

  }

  // A function that will tell if a track is located a distance 'd' from the edge of the active volume of the TPC.
  // Inputs: pt - the object of type '::geoalgo::Vector' that is being checked for being within a distance 'd' of a detector boundary.
  //         d  - the distance from the detector boundary for which the input point is being checked.
  bool GeneralFlashMatchAlgo::isNearActiveVolumeEdge( ::geoalgo::Vector pt, double d ) {

    // Return 'true' if the input point is near one of the edges of the detector.
    if ( pt[0] <= (0.0 + d ) || pt[0] >= (256.35 - d ) || pt[1] <= (-116.5 + d) || pt[1] >= (116.5 - d) || pt[2] <= (0.0 + d) || pt[2] >= (1036.8 - d) ) {
      return true;
    }

    // Return false otherwise.
    return false;
  }
  
  // A function that will find the center and range of a flash.
  // Inputs: flash: This is hte larlite flash object that is being considered.
  //         zmean: This is the mean z position of the flash.
  //         zrange: This is the range in the z variable that we are considering.
  void GeneralFlashMatchAlgo::GetFlashCenterAndRange( const larlite::opflash& flash, float& zmean, std::vector<float>& zrange ) {

    zmean = 0.;
    zrange.resize(2,0.0);
    float qtot = 0.0;
    for (int ipmt=0; ipmt<32; ipmt++) {
      zmean += flash.PE(ipmt)*m_pmtweights.pmtpos[ipmt][2];
      qtot  += flash.PE( ipmt );
    }

    if (qtot>0)
      zmean /= qtot;

    float mindist = 1000;
    float maxdist = 0;

    for (int ipmt=0; ipmt<32; ipmt++) {
      if (flash.PE(ipmt)>m_config.pmtflash_thresh) {
        float dist = m_pmtweights.pmtpos[ipmt][2]-zmean;
        if (dist<0 || mindist>dist)
          mindist = dist;
        else if ( dist>0 && maxdist<dist )
          maxdist = dist;
      }
    }

    if (qtot<=0 )  {
      // no information, so set wide range                                                                                                                                                                  
      zrange[0] = 0;
      zrange[1] = 3455;
    }
    else {
      zrange[0] = (zmean+mindist-10.0)/0.3; // conversion from cm -> wire                                                                                                                                   
      zrange[1] = (zmean+maxdist+10.0)/0.3; // conversion from cm -> wire                                                                                                                                   
    }

    // bounds check                                                                                                                                                                                         
    if ( zrange[0]<0 ) zrange[0] = 0;
    if ( zrange[1]>=3455 ) zrange[1] = 3455;

  }

  // A function that will generate an unfitted flash hypothesis for the track.
  // Inputs: qcluster: This is the qcluster that will be used to generate the flash hypothesis.
  flashana::Flash_t GeneralFlashMatchAlgo::GenerateUnfittedFlashHypothesis( const flashana::QCluster_t& qcluster ) {
    const flashana::QLLMatch* matchalgo = (flashana::QLLMatch*)m_flash_matcher.GetAlgo( flashana::kFlashMatch );
    flashana::Flash_t unfitted_hypothesis = matchalgo->GetEstimate( qcluster );
    return unfitted_hypothesis;
  }
  
  // A function that will make an opflash from a data flash.
  // Input: Flash: This is an object of type 'Flash_t' (a data flash).
  larlite::opflash GeneralFlashMatchAlgo::MakeOpFlashFromFlash(const flashana::Flash_t& Flash) {

    const larutil::Geometry* geo = ::larutil::Geometry::GetME();

    std::vector<double> PEperOpDet( Flash.pe_v.size(), 0.0);
    for (unsigned int femch = 0; femch < Flash.pe_v.size(); femch++) {
      unsigned int opdet = geo->OpDetFromOpChannel(femch);
      PEperOpDet[femch] = Flash.pe_v[opdet];
    }

    larlite::opflash flash( Flash.time, 0, 0, 0, PEperOpDet, false, false, 1.0, Flash.y, Flash.y_err, Flash.z, Flash.z_err );

    return flash;
  }

  // A function that will calculate the 2D gaussian value for an opflash object.
  // Inputs: 'flash' - The input 'opflash' object to be fit to a Gaussian 2D value.
  larlite::opflash GeneralFlashMatchAlgo::GetGaus2DPrediction( const larlite::opflash& flash ) {
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();
    float tot_weight = 0;

    // Get Charge Weighted Mean
    float mean[2] = {0.0};
    std::vector<double> xyz(3,0.0);
    for (int iopdet=0; iopdet<32; iopdet++ ) {
      larutil::Geometry::GetME()->GetOpDetPosition( iopdet, xyz );
      mean[0] += xyz[1]*flash.PE(iopdet);
      mean[1] += xyz[2]*flash.PE(iopdet);
      tot_weight += flash.PE(iopdet);
    }
    for (int i=0; i<2; i++)
      mean[i] /= tot_weight;

    /// Get Charge-weighted covariance Matrix
    float cov[2][2] = {0};
    for (int iopdet=0; iopdet<32; iopdet++) {
      larutil::Geometry::GetME()->GetOpDetPosition( iopdet, xyz );
      float dx[2];
      dx[0] = xyz[1] - mean[0];
      dx[1] = xyz[2] - mean[1];
      for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++ ) {
          cov[i][j] += dx[i]*dx[j]*flash.PE(iopdet)/tot_weight;
        }
      }
    }

    // Get Inverse Covariance Matrix
    float det = cov[0][0]*cov[1][1] - cov[1][0]*cov[0][1];
    float invcov[2][2] = { { cov[1][1]/det, -cov[0][1]/det},
                           {-cov[1][0]/det,  cov[0][0]/det } };
    float normfactor = 1.0/sqrt(2*3.141159*det);


    // build guas ll flash hypothesis using cov matrix
    std::vector< double > PEperOpDet(32,0);
    float normw = 0.;
    for ( int iopdet=0; iopdet<32; iopdet++) {
      larutil::Geometry::GetME()->GetOpDetPosition( iopdet, xyz );
      float mahadist = 0.;
      float yz[2] = { (float)(xyz[1]), (float)(xyz[2]) };
      for (int i=0; i<2; i++ ) {
        for (int j=0; j<2; j++) {
          mahadist += (yz[i]-mean[i])*invcov[i][j]*(yz[j]-mean[j]);
        }
      }
      PEperOpDet[iopdet] = normfactor*exp(-0.5*mahadist);
      normw += PEperOpDet[iopdet];
    }

    // normalize back to totalweight                                            
    if ( normw>0 ) {
      for (int iopdet=0; iopdet<32; iopdet++) {
        PEperOpDet[iopdet] *= tot_weight/normw;
      }
    }

    // make flash                                                               
    larlite::opflash gaus2d_flash( flash.Time(), flash.TimeWidth(), flash.AbsTi\
me(), flash.Frame(), PEperOpDet );
    return gaus2d_flash;
  }


  // A function that will make a data flash (an object of type 'Flash_t') from an 'opflash'.
  // Inputs: flash: This is an object of type 'opflash' that has to be converted into a data flash.
  flashana::Flash_t GeneralFlashMatchAlgo::MakeDataFlash( const larlite::opflash& flash ) {

    // get geometry info                                                                                                                                                                                    
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();

    // provide mean and basic width                                                                                                                                                                         
    float wire_mean;
    std::vector<float> wire_range;
    GetFlashCenterAndRange( flash, wire_mean, wire_range );
    if ( m_verbosity>0 )
      std::cout << "Flash: pe=" << flash.TotalPE() << " wire_mean=" << wire_mean << " wire_range=[" << (int)wire_range[0] << "," << (int)wire_range[1] << "]" << std::endl;

    // also load flash info into flash manager                                                                                                                                                              
    flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i) / m_config.gain_correction[i];
      f.pe_err_v[opdet] = sqrt(flash.PE(i) / m_config.gain_correction[i]);
    }
    f.time = flash.Time();
    f.idx = 0;

    return f;
  }


  // Write a function from the 'TaggerFlashMatchAlgo' that returns an entire vector of data flashes.
  std::vector<flashana::Flash_t> GeneralFlashMatchAlgo::MakeDataFlashes( std::vector<larlite::opflash> opflash_v ) {
    std::vector<flashana::Flash_t> flashes_v;
    for ( auto const& opflash : opflash_v ) {
      flashana::Flash_t flash = MakeDataFlash( opflash );
      flashes_v.emplace_back( std::move(flash) );
    }
    return flashes_v;
  }
  
  // A function that will generate a flash hypothesis for a larlite track.  This will follow the logic at the beginning of the 'GeneralFlashMatchAlgo::InTimeFlashComparison' function.
  // Inputs: input_track: This is the input larlite track needed to generate the flash hypothesis.  It will be converted to a qcluster, which will be converted a flash_t object, which will be
  // converted into an opflash object.
  // *** We will not be using the gaus2d object here right now. ****
  //         taggerflashmatch_pset: This is the input parameter set used to instatiate the GeneralFlashMatchAlgo object.
  larlite::opflash GeneralFlashMatchAlgo::make_flash_hypothesis(const larlite::track input_track) {

    // Declare a qcluster for use in the 'ExpandQClusterStartingWithLarliteTrack' function.
    flashana::QCluster_t qcluster;

    // Follow the logic of the 'GeneralFlashMatchAlgo::InTimeFlashComparison' function to generate a flash hypothesis.
    ExpandQClusterStartingWithLarliteTrack(qcluster, input_track, 10000., true, true); 
    flashana::Flash_t flash_hypo  = GenerateUnfittedFlashHypothesis( qcluster );
    larlite::opflash opflash_hypo = MakeOpFlashFromFlash( flash_hypo );
    larlite::opflash gaus2d_hypo  = GetGaus2DPrediction( opflash_hypo );

    // What to return depends on the values set in the 'Config' class.
    if (!m_config.use_gaus2d)
      return opflash_hypo;

    else
      return gaus2d_hypo;

  }

  
  // A function that will make flash hypotheses from an entire vector of larlite tracks.
  std::vector<larlite::opflash> GeneralFlashMatchAlgo::make_flash_hypothesis_vector( const std::vector<larlite::track>& input_track_v ) {
    std::vector<larlite::opflash> opflash_hypo_v;
    for ( auto& input_track : input_track_v ) {
      larlite::opflash opflash_hypo = make_flash_hypothesis( input_track );
      opflash_hypo_v.emplace_back( std::move(opflash_hypo) );
    }
    return opflash_hypo_v;
  }
  
// A function that will take in reconstructed track and flash info and return a chi2 fit between them, using the 'GeneralFlashMatchAlgo::InTimeFlashComparison' class.
// Inputs: opflash_hypo: This is the flash hypothesis that is used in compare to the data flash.
//         qcluster: The qcluster that corresponds to the opflash hypothesis.
//         data_opflash: This is the 'data_flash' that you are matching to the 'flash_hypothesis' opflash object for the track.
//         flash_prod_idx: This is a binary value for the flash producer which determines if we correct the waveform. '0' for the 'simpleFlashBeam' (no correction) or '1' for 'simpleFlashCosmic' (with correction).
  float GeneralFlashMatchAlgo::generate_chi2_in_track_flash_comparison(const flashana::QCluster_t qcluster, const larlite::opflash data_opflash, float& totpe_data_flash, float& totpe_hypo_flash, int flash_prod_idx) {
    
    // Convert 'data_opflash' into type 'Flash_t' (a data flash).
    flashana::Flash_t data_flasht = MakeDataFlash( data_opflash );

    // Convert 'qcluster' into an unfitted flash hypothesis.
    flashana::Flash_t flash_hypo = GenerateUnfittedFlashHypothesis( qcluster );

    // Convert the 'flash_hypo' to a data flash and then a Gaussian 2D object to have that below.
    larlite::opflash opflash_hypo = MakeOpFlashFromFlash( flash_hypo );
    larlite::opflash gaus2d_hypo  = GetGaus2DPrediction( opflash_hypo );
    
    // Use the geometry package at this point in the algorithm.
    const larutil::Geometry* geo = ::larutil::Geometry::GetME();

    // Declare the variables for finding the chi2 value.
    float chi2       = 0.;
    totpe_data_flash = 0.;
    totpe_hypo_flash = 0.;
    std::vector< float > expectation(32,0);

    // Declare a vector for the number of valid PMTs.
    int nvalid_pmts = 0;

    // This begins the part of the code that finds the chi2 value.
    for (size_t i=0; i<data_flasht.pe_v.size(); i++) {
      float observed = data_flasht.pe_v.at(i);
      
      // Apply these changes IF the flash is outside of the unbiased PMT (the beam) window.
      if ( flash_prod_idx == 1 ) {

	if ( flash_hypo.pe_v.at( i ) < 60.0 ) 
	  flash_hypo.pe_v.at( i ) *= 0.424;

	else if ( flash_hypo.pe_v.at( i ) > 60.0 )
	  flash_hypo.pe_v.at( i) *= 0.354;

      }
      
      // Declare 'expected' down here after the cosmic discriminator correction has been applied.
      float expected    = flash_hypo.pe_v.at(i);

      totpe_data_flash += observed;
      totpe_hypo_flash += expected;

      if ( observed>0 && expected==0 ) {
	if ( !m_config.use_gaus2d )
	  expected = 1.0e-3;
	else
	  expected = gaus2d_hypo.PE( geo->OpDetFromOpChannel(i) )*m_config.fudge_factor;
      }
      expectation[i] = expected;

      float dpe = 0.;
      if ( observed>0 ) {
	dpe   = (expected-observed) + observed*( log( observed ) - log( expected ) );
	chi2 += 2.0*dpe;
	++nvalid_pmts;
      }
    }
    chi2 /= nvalid_pmts;

    std::cout << "  [observed] ";
    for (int ich=0; ich<32; ich++ ) {
      std::cout << std::setw(5) << (int)data_flasht.pe_v.at(  geo->OpDetFromOpChannel(ich) );
    }
    std::cout << " TOT=" << totpe_data_flash << " CHI2=" << chi2 << std::endl;

    std::cout << "  [expected] ";
    for ( int ich=0; ich<32; ich++ ) {
      //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at(  geo->OpDetFromOpChannel(ich) )*m_config.fudge_factor);
      //std::cout << std::setw(5) << (int)(opflash_hypo.pe_v.at( ich )*m_config.fudge_factor);
      std::cout << std::setw(5) << (int)expectation[ich];
    }
    std::cout << " TOT=" << totpe_hypo_flash << " BestCHI2=" << chi2 << std::endl;

    // Return the chi2 value.
    return chi2;
        
  }
 
}

