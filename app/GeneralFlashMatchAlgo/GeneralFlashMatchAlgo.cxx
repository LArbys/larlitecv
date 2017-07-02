// (All of the functions for flash matching from 'ContainedROI/TaggerFlashMatchAlgo.h' are included in "GeneralFlashMatchAlgo.h".
#include "GeneralFlashMatchAlgo.h"

namespace larlitecv {

  // A function that will find the flash indices that correspond to flashes that do not correspond to anode/cathode piercing endpoints.
  // Inputs: opflash_v: all of the opflash objects from the event.
  //         anode_flash_idx_v: the indices of the flashes that correspond to boundary points at the anode.
  //         cathode_flash_idx_v: the indices of the flashes that correspond to boundary points at the cathode.
  std::vector < int > GeneralFlashMatchAlgo::non_anodecathode_flash_idx( const std::vector < larlite::opflash* >& opflash_v, std::vector < int > anode_flash_idx_v, std::vector < int > cathode_flash_idx_v ) {

    // Declare a vector for the flashes that are neither anode or cathode crossing.
    std::vector < int > non_ac_flash_idx;

    // Loop through the opflash information from the event in order to separate the anode/cathode crossing flashes from the others.
    for (size_t non_ac_i = 0; non_ac_iter < opflash_v.size(); non_ac_iter++) {

      bool belongs_to_ac_flash = false;

      // Loop through the anode flashes.
      for (size_t anode_flsh_i = 0; anode_flsh_i < anode_flash_idx_v.size(); anode_flsh_i++) {

	// Compare the two flash indices.
	if (non_ac_i == anode_flash_idx_v[anode_flsh_i]) belongs_to_ac_flash = true;

      }

      // Loop through the cathode flashes.
      for (size_t cathode_flsh_i = 0; cathode_flsh_i < cathode_flash_idx_v.size(); cathode_flsh_i++) {

        // Compare the two flash indices.
	if (non_ac_i == cathode_flash_idx_v[cathode_flsh_i]) belongs_to_ac_flash = true;

      }

      // If it was not the same index value as any of the anode/cathode flash indices, then you can append it onto the 'non_ac_flash_idx' vector.
      if (!belongs_to_ac_flash) non_ac_flash_idx.push_back(non_ac_i);

    }

    // Return the vector of non anode-piercing/cathode-piercing track endpoints.
    return non_ac_flash_idx;

  }

  // A function that will generate a list of opflashes that correspond to the list of indices that are put in.
  // Inputs: full_opflash_v: This is the full list of opflashes from which we are trying to filter opflashes of the same denomination, i.e. all anode-piercing flashes, all cathode-piercing flashes, all non-side-piercing flashes, etc.
  //         idx_v: This is the list of indices of the flashes that you are interested in.  You want to place the flashes located in the 'full_opflash_v' vector at each index in this list in the output list, which contains the opflashes that we are interested in.
  std::vector <opflash*> GeneralFlashMatchAlgo::generate_single_denomination_flash_list(const std::vector< larlite::opflash > full_opflash_v, std::vector < int > idx_v) {

    // Declare a list of flashes that will be the output, which is the list of the certain denomination of flashes.
    std::vector < larlite::opflash > single_denomination_flash_list;
    single_denomination_flash_list.clear();
    
    for (size_t iflash = 0; iflash < idx_v.size(); iflash++) {

      single_denomination_flash_list.push_back(full_opflash_v.at(idx.at(iflash)));

    }

    // Return this list of the opflash information of one denomination.
    return single_denomination_flash_list;

  }

  // A function that will generate tracks from the 'BMTrackCluster3D' objects that are formed from each pass of the tagger.
  // Inputs: trackcluster3d_v: A vector of 'BMTrackCluster3D' objects that were formed from a pass of the tagger.
  std::vector < larlite::track >  GeneralFlashMatchAlgo::generate_tracks_between_passes(std::vector< BMTrackCluster3D > trackcluster3d_v) {

    // We will use the 'makeTrack()' functionality of these 'BMTrackCluster3D' objects to turn them into a vector of 'larlite::track' objects.
    std::vector < larlite::track > output_tracks;

    // Loop through the list of 'BMTrackCluster3D' objects and turn them into a vector 'larlite::track' objects.
    for (auto const& trackcluster3d : trackcluster3d_v) {

      output_tracks.push_back(trackcluster3d.makeTrack());

    }

    // Return this list of tracks.
    return output_tracks;

  }

  // A function that will generate a flash hypothesis for a larlite track.  This will follow the logic at the beginning of the 'TaggerFlashMatchAlgo::InTimeFlashComparison' function.
  // Inputs: input_track: This is the input larlite track needed to generate the flash hypothesis.  It will be converted to a qcluster, which will be converted a flash_t object, which will be
  // converted into an opflash object.
  // *** We will not be using the gaus2d object here right now. ****
  //         taggerflashmatch_pset: This is the input parameter set used to instatiate the TaggerFlashMatchAlgo object.
  larlite::opflash make_flash_hypothesis(const larlite::track input_track, larcv::PSet taggerflashmatch_pset) {

    larlitecv::TaggerFlashMatchAlgoConfig taggerflashmatch_cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( taggerflashmatch_pset );
    larlitecv::TaggerFlashMatchAlgo taggerflashmatch( taggerflashmatch_cfg );

    // Follow the logic of the 'TaggerFlashMatchAlgo::InTimeFlashComparison' function to generate a flash hypothesis.
    flashana::QCluster_t qcluster = taggerflashmatch.GenerateQCluster(input_track);
    flashana::Flash_t flash_hypo  = taggerflashmatch.GenerateUnfittedFlashHypothesis( qcluster );
    larlite::opflash opflash_hypo = taggerflashmatch.MakeOpFlashFromFlash( flash_hypo );

    // I will not generate the gaus2d object within this function, which is done in the 'TaggerFlashMatchAlgo' class of 'ContainedROI'.
    return opflash_hypo;

  }
    
// A function that will take in reconstructed track and flash info and return a chi2 fit between them, using the 'TaggerFlashMatchAlgo::InTimeFlashComparison' class.
// Inputs: input_track: This is the input track needed for the flash-matching.
//         data_opflash: This is the 'data_flash' that you are matching to the 'flash_hypothesis' opflash object for the track.
//         taggerflashmatch_pset: This is the parameter set containing all of the parameters to set up the flash matching for the ContainedROI stage.
  float GeneralFlashMatchAlgo::generate_chi2_in_track_flash_comparison(const larlite::track input_track, const larlite::opflash data_opflash, larcv::PSet taggerflashmatch_pset) {
    
    // Using the given parameter set for the ContainedROI class, declare the objects needed to turn the input types (qcluster, opflash) into those that are needed by the 'TaggerFlashMatchAlgo::InTimeFlashComparison' function (opflash (the hypothesis), std::vector <Flash_t>).
    larlitecv::TaggerFlashMatchAlgoConfig taggerflashmatch_cfg = larlitecv::TaggerFlashMatchAlgoConfig::MakeTaggerFlashMatchAlgoConfigFromPSet( taggerflashmatch_pset );
    larlitecv::TaggerFlashMatchAlgo taggerflashmatch( taggerflashmatch_cfg );
    
    // Convert 'input_track' into type 'qcluster'.
    flashana::QCluster_t qcluster = taggerflashmatch.GenerateQCluster(input_track);

    // Convert 'data_opflash' into type 'Flash_t' (a data flash).
    flashana::Flash_t data_flasht = MakeDataFlash(data_opflash);

    // Convert this into a vector, which is what the 'TaggerFlashMatchAlgo::InTimeFlashComparison' takes as input.
    std::vector < flashana::Flash_t > data_flasht_v;
    data_flasht_v.resize(1, 0.0);

    // Initialize the vector by appending the 'data_flasht' object onto the end.
    data_flasht_v[0] = data_flasht;

    // Declare the two inputs needed for the
    float tot_pe_hypo = 0.0;
    float totpe_data  = 0.0;

    // Use these inputs in the 'TaggerFlashMatchAlgo::InTimeFlashComparison' function to obtain a chi2 fit.
    float chi2 = taggerflashmatch.InTimeFlashComparison(data_flasht_v, qcluster, totpe_data, tot_pe_hypo);

    // Return the chi2 value.
    return chi2;
        
  }

}
