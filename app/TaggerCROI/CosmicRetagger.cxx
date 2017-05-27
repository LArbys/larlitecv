#include "CosmicRetagger.h"

// larcv
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"

// larlitecv
#include "GapChs/EmptyChannelAlgo.h"
#include "UnipolarHack/UnipolarHackAlgo.h"
#include "TaggerCROIAlgo.h"

namespace larlitecv {

  CosmicRetagger::CosmicRetagger( std::string retagger_cfg, std::string taggerout_larcv_file, std::string taggerout_larlite_file ) {
    setConfigFile( retagger_cfg );

    // we have a custom IOManager configuration. It reads every tree as we don't know what trees were written by the version of the tagger un.
    //std::string larcv_io_cfg=" Verbosity: 2 IOMode: 0 InputFiles: [] ReadOnlyDataType: [] ReadOnlyDataName: []";
    //std::string larlite_io_cfg=" Verbosity: 2 IOMode: 0 ReadOnlyProducers: [] ReadOnlyDataTypes: []";

    std::string larcv_io_cfg=" Verbosity: 2 IOMode: 0";
    std::string larlite_io_cfg=" Verbosity: 2 IOMode: 0";
    
    larcv::PSet larcv_pset   = larcv::PSet( "IOManager", larcv_io_cfg );
    larcv::PSet larlite_pset = larcv::PSet( "StorageManager", larlite_io_cfg );

    // input
    m_dataco_input.configure( larcv_pset, larlite_pset );
    m_dataco_input.add_inputfile( taggerout_larcv_file, "larcv" );
    m_dataco_input.add_inputfile( taggerout_larlite_file, "larlite" );
    m_dataco_input.initialize();

    // output
    m_dataco_output.configure( m_cfg_file, "StorageManagerOut", "IOManagerOut", "TaggerCROI" );
    m_dataco_output.initialize();

    m_nentries = m_dataco_input.get_nentries("larcv");
    std::cout << "CosmicRetagger Entries: larcv=" << m_nentries << " larlite=" << m_dataco_input.get_nentries("larlite") << std::endl;

  }

  bool CosmicRetagger::processInputImages() {
    // we look into the file and check what's inside
    if ( !m_state.configured ) {
      std::cout << "Invalid state to run processInputImages: " << printState() << std::endl;
      return false;
    }
    
    // set the state: we reset all downstream states
    m_state.input_ready = m_state.thrumu_run = m_state.stopmu_run = m_state.untagged_run = m_state.croi_run = false;
    
    // -----------------------------------------------------------------------------
    // Load the data
    int run, subrun, event;
    try {
      m_dataco_input.goto_entry( m_entry, "larcv" );

      // set the id's
      m_dataco_input.get_id(run,subrun,event);
      m_dataco_output.set_id(run,subrun,event);
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading the input data: " << e.what() << std::endl;
      return false;
    }
      
    std::cout << "=====================================================" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << " Entry " << m_entry << " : " << run << " " << subrun << " " << event << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // get images (from larcv)
    m_input_data.clear();
    larcv::EventImage2D* event_imgs    = (larcv::EventImage2D*)m_dataco_input.get_larcv_data( larcv::kProductImage2D, "modimg" );
    if ( event_imgs->Image2DArray().size()==0) {
      if ( !m_skip_empty_events )
        throw std::runtime_error("Number of images=0. LArbys.");
      else {
        std::cout << "Skipping Empty Events." << std::endl;
	return false;
      }
    }
    else if ( event_imgs->Image2DArray().size()!=3 ) {
      throw std::runtime_error("Number of Images!=3. Weird.");
    }
    event_imgs->Move( m_input_data.img_v );
    const std::vector<larcv::Image2D>& imgs_v = m_input_data.img_v;
    
    // Reload ThruMu
    // reload ThruMu payload data: need track points and pixels. and tagged image.
    
    larcv::EventPixel2D*  ev_thrumu_pixels = (larcv::EventPixel2D*)m_dataco_input.get_larcv_data( larcv::kProductPixel2D, "thrumupixels" );
    larlite::event_track* ev_thrumu_tracks = (larlite::event_track*)m_dataco_input.get_larlite_data( larlite::data::kTrack, "thrumu3d" );
    if ( ev_thrumu_tracks==NULL )
      throw std::runtime_error("Could not load thrumu tracks from larlite file.");
    if ( ev_thrumu_tracks->size()!=ev_thrumu_pixels->Pixel2DClusterArray((larcv::PlaneID_t)0).size() ) {
      std::cout << "number of tracks (" << ev_thrumu_tracks->size() << ")"
		<< " does not match number of pixel clusters (" << ev_thrumu_pixels->Pixel2DClusterArray((larcv::PlaneID_t)0).size() << ")" << std::endl;
      throw std::runtime_error("Not good.");
    }

    m_thrumu_data.clear();
    
    for (size_t itrack=0; itrack<ev_thrumu_tracks->size(); itrack++) {
      m_thrumu_data.track_v.push_back( ev_thrumu_tracks->at(itrack) );
      std::vector<larcv::Pixel2DCluster> pixelcluster_v;
      for ( size_t p=0; p<imgs_v.size(); p++ ) {
        const larcv::Pixel2DCluster& pixels = ev_thrumu_pixels->Pixel2DClusterArray((larcv::PlaneID_t)p).at(itrack);
        pixelcluster_v.push_back( pixels );
      }
      m_thrumu_data.pixelcluster_v.emplace_back( std::move(pixelcluster_v) );
    }

    m_thrumu_data.tagged_v.clear();
    for ( size_t p=0; p<imgs_v.size(); p++ ) {
      larcv::Image2D thrumu_img( imgs_v.at(p).meta() );
      thrumu_img.paint(0);
      for ( auto const& pixelcluster_v : m_thrumu_data.pixelcluster_v ) {
        for ( auto const& pixel : pixelcluster_v.at(p) ) {
          for (int dr=-5; dr<=5; dr++) {
            int row = dr + (int)pixel.Y();
            if ( row<0 || row>=(int)imgs_v.at(p).meta().rows() ) continue;
            for (int dc=-5; dc<5; dc++) {
              int col=dc+(int)pixel.X();
              if ( col<0 || col>=(int)imgs_v.at(p).meta().cols() ) continue;
              thrumu_img.set_pixel(row,col,255);
            }//end of col loop
          }//end of row loop
        }//end of pixel loop
      }//end of cluster loop
      m_thrumu_data.tagged_v.emplace_back( std::move(thrumu_img) );
    }//end of plane loop
    
    return true;
  }
  
}
