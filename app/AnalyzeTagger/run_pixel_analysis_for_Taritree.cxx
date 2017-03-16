#include "Base/DataCoordinator.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// CRB: For the information for the pixels in each event.
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Pixel2DCluster.h"
#include "DataFormat/Pixel2D.h"

// larlite
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "TFile.h"
#include "TTree.h"

float dwall( const std::vector<float>& pos, int& boundary_type ) {

  float dx1 = fabs(pos[0]);
  float dx2 = fabs(255-pos[0]);
  float dy1 = fabs(116.0-pos[1]);
  float dy2 = fabs(-116.0-pos[1]);
  float dz1 = fabs(pos[2]);
  float dz2 = fabs(1036.0-pos[2]);

  float dwall = 1.0e9;

  if ( dy1<dwall ) {
  	dwall = dy1;
  	boundary_type = 0;
  }
  if ( dy2<dwall ) {
  	dwall = dy2;
  	boundary_type = 1;
  }
  if ( dz1<dwall ) {
  	dwall = dz1;
  	boundary_type = 2;
  }
  if ( dz2<dwall ) {
  	dwall = dz2;
  	boundary_type = 3;
  }
  if ( dx1<dwall ) {
  	dwall = dx1;
  	boundary_type = 4;
  }
  if ( dx2<dwall ) {
  	dwall = dx2;
  	boundary_type = 5;
  }

  return dwall;

}


int main( int nargs, char** argv ) {

  // Declaring the function 'isNeutrinoPixel' right here.
  bool isNeutrinoPixel(const larcv::Image2D& segmented_image, float pix_x, float pix_y);

  // run pixel analysis. set the config parameters first and then you can place them in the files below.
  // configuration parameters                                                                                                                                           
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet pixana_cfg = cfg.get<larcv::PSet>("PixelAnalysis");
  float fthreshold   = pixana_cfg.get<float>("PixelThreshold");                                                                                                                  
  int verbosity      = pixana_cfg.get<int>("Verbosity",0);

  larlitecv::DataCoordinator dataco_source;
  dataco_source.set_filelist( "taggeranalysis_original_img_input_file_list_larcv.txt", "larcv" ); // segment image/original image
  dataco_source.set_filelist( "taggeranalysis_original_img_input_file_list_larlite.txt", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
  dataco_thrumu.set_filelist( "taggeranalysis_thrumu_input_file_list_larcv.txt", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.set_filelist( "taggeranalysis_thrumu_input_file_list_larlite.txt", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.set_filelist( "taggeranalysis_stopmu_input_file_list_larcv.txt" , "larcv" ); //stopmu-tagger output

  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_thrumu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_stopmu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

  // setup output

  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
  TTree* tree = new TTree("pixana", "Pixel-level analysis");
  int run, subrun, event;

  // Declare variables for the quantities that I am interested in plotting.  These are for the quantities including ALL of the planes.
  float _frac_all_pixels_tagged_more_than_once;
  float _frac_cosmic_pixels_tagged_once_or_more;
  float _frac_cosmic_pixels_tagged_more_than_once;
  float _frac_cosmic_pixels_not_tagged;
  float _frac_neutrino_pixels_not_tagged;

  float _all_planes_all_pixels_tree;
  float _all_planes_cosmic_pixels_tree;
  float _all_planes_neutrino_pixels_tree;

  // Declare the same information but specific to each plane.
  // u
  float _u_plane_frac_all_pixels_tagged_more_than_once;
  float _u_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _u_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _u_plane_frac_cosmic_pixels_not_tagged;
  float _u_plane_frac_neutrino_pixels_not_tagged;

  // Include the number of each type of pixel used on the plane.
  float _u_plane_all_pixels;
  float _u_plane_cosmic_pixels;
  float _u_plane_neutrino_pixels;

  // v
  float _v_plane_frac_all_pixels_tagged_more_than_once;
  float _v_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _v_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _v_plane_frac_cosmic_pixels_not_tagged;
  float _v_plane_frac_neutrino_pixels_not_tagged;

  float _v_plane_all_pixels;
  float _v_plane_cosmic_pixels;
  float _v_plane_neutrino_pixels;

  // y
  float _y_plane_frac_all_pixels_tagged_more_than_once;
  float _y_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _y_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _y_plane_frac_cosmic_pixels_not_tagged;
  float _y_plane_frac_neutrino_pixels_not_tagged;

  float _y_plane_all_pixels;
  float _y_plane_cosmic_pixels;
  float _y_plane_neutrino_pixels;

  // Declare a value for the current pixel's value.
  float _current_all_pixel_value;
  float _current_neutrino_pixel_value;
  float _current_cosmic_pixel_value;

  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");

  // Create branches for the quantities that I am interested in plotting.
  tree->Branch("_frac_all_pixels_tagged_more_than_once", &_frac_all_pixels_tagged_more_than_once, "frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_frac_cosmic_pixels_tagged_once_or_more", &_frac_cosmic_pixels_tagged_once_or_more, "frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_frac_cosmic_pixels_tagged_more_than_once", &_frac_cosmic_pixels_tagged_more_than_once, "frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_frac_cosmic_pixels_not_tagged", &_frac_cosmic_pixels_not_tagged, "frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_frac_neutrino_pixels_not_tagged", &_frac_neutrino_pixels_not_tagged, "frac_neutrino_pixels_not_tagged/F");

  tree->Branch("_all_planes_all_pixels", &_all_planes_all_pixels_tree, "all_planes_all_pixels/F");
  tree->Branch("_all_planes_cosmic_pixels", &_all_planes_cosmic_pixels_tree, "all_planes_cosmic_pixels/F");
  tree->Branch("_all_planes_neutrino_pixels", &_all_planes_neutrino_pixels_tree, "all_planes_neutrino_pixels/F");

  // Create branches for the same quantities but for each of the wire planes.
  // u
  tree->Branch("_u_plane_frac_all_pixels_tagged_more_than_once", &_u_plane_frac_all_pixels_tagged_more_than_once, "u_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_tagged_once_or_more", &_u_plane_frac_cosmic_pixels_tagged_once_or_more, "u_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_tagged_more_than_once", &_u_plane_frac_cosmic_pixels_tagged_more_than_once, "u_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_not_tagged", &_u_plane_frac_cosmic_pixels_not_tagged, "u_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_u_plane_frac_neutrino_pixels_not_tagged", &_u_plane_frac_neutrino_pixels_not_tagged, "u_plane_frac_neutrino_pixels_not_tagged/F");

  tree->Branch("_u_plane_all_pixels", &_u_plane_all_pixels, "u_plane_all_pixels/F");
  tree->Branch("_u_plane_cosmic_pixels", &_u_plane_cosmic_pixels,"u_plane_cosmic_pixels/F");
  tree->Branch("_u_plane_neutrino_pixels", &_u_plane_neutrino_pixels, "u_plane_neutrino_pixels/F");

  // v
  tree->Branch("_v_plane_frac_all_pixels_tagged_more_than_once", &_v_plane_frac_all_pixels_tagged_more_than_once, "v_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_tagged_once_or_more", &_v_plane_frac_cosmic_pixels_tagged_once_or_more, "v_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_tagged_more_than_once", &_v_plane_frac_cosmic_pixels_tagged_more_than_once, "v_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_not_tagged", &_v_plane_frac_cosmic_pixels_not_tagged, "v_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_v_plane_frac_neutrino_pixels_not_tagged", &_v_plane_frac_neutrino_pixels_not_tagged, "v_plane_frac_neutrino_pixels_not_tagged/F");
  
  tree->Branch("_v_plane_all_pixels", &_v_plane_all_pixels, "v_plane_all_pixels/F");
  tree->Branch("_v_plane_cosmic_pixels", &_v_plane_cosmic_pixels,"v_plane_cosmic_pixels/F");
  tree->Branch("_v_plane_neutrino_pixels", &_v_plane_neutrino_pixels, "v_plane_neutrino_pixels/F");

  // y
  tree->Branch("_y_plane_frac_all_pixels_tagged_more_than_once", &_y_plane_frac_all_pixels_tagged_more_than_once, "y_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_tagged_once_or_more", &_y_plane_frac_cosmic_pixels_tagged_once_or_more, "y_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_tagged_more_than_once", &_y_plane_frac_cosmic_pixels_tagged_more_than_once, "y_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_not_tagged", &_y_plane_frac_cosmic_pixels_not_tagged, "y_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_y_plane_frac_neutrino_pixels_not_tagged", &_y_plane_frac_neutrino_pixels_not_tagged, "y_plane_frac_neutrino_pixels_not_tagged/F");

  tree->Branch("_y_plane_all_pixels", &_y_plane_all_pixels, "y_plane_all_pixels/F");
  tree->Branch("_y_plane_cosmic_pixels", &_y_plane_cosmic_pixels,"y_plane_cosmic_pixels/F");
  tree->Branch("_y_plane_neutrino_pixels", &_y_plane_neutrino_pixels, "y_plane_neutrino_pixels/F");


  int nentries = dataco_stopmu.get_nentries("larcv");

  dataco_source.goto_entry(0,"larcv");
  dataco_thrumu.goto_entry(0,"larcv");
  dataco_stopmu.goto_entry(0,"larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    // 'dataco_thrumu' entries.
    dataco_thrumu.goto_entry(ientry,"larcv");
    
    dataco_thrumu.get_id(run,subrun,event);

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << std::endl;
      std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    dataco_thrumu.goto_event(run,subrun,event,"larcv");
    dataco_source.goto_event(run,subrun,event,"larcv");

    // Add a line for 'stopmu' (not sure why it isn't already there).
    dataco_stopmu.goto_event(run,subrun,event,"larcv");

    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_segs   = (larcv::EventImage2D*)dataco_source.get_larcv_data(larcv::kProductImage2D,"seg_comb_tpc"); // Previously 'segment' in the second entry.
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"modimgs"); // Ask Taritree about this, could be wrong file.
    larcv::EventImage2D* ev_badch  = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"gapchs");  // Change the information in here
    larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"marked3d");
    larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D,"stopmu");
    
    // Do the analogous thing for the track information here.
    // This is 'truth' information, so this should also be from the 'data_source' file.
    larcv::EventPixel2D* ev_thrumu_pixel_clusters = (larcv::EventPixel2D*)dataco_thrumu.get_larcv_data(larcv::kProductPixel2D, "thrumu2d"); 
    larcv::EventPixel2D* ev_stopmu_pixel_clusters = (larcv::EventPixel2D*)dataco_stopmu.get_larcv_data(larcv::kProductPixel2D, "stopmupixels");

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>& segs_v   = ev_segs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
    const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();

    // Declare variables for the total number of both types of pixels, and pixels in general, on each of the planes.
    // These data types will have to be of type 'float' so that they can be used in the fractional values that I am looking at.

    // Declare variables for the 'truth' number of pixels.
    float _all_planes_cosmic_pixels          = 0.;
    float _all_planes_neutrino_pixels        = 0.;
    float _all_planes_all_pixels             = 0.;

    // Declare variables for the 'tagged' number of pixels.
    float _all_planes_cosmic_pixels_tagged   = 0.;
    float _all_planes_neutrino_pixels_tagged = 0.;
    float _all_planes_all_pixels_tagged      = 0.;

    // Declare variables for the number of pixels (of each type) tagged once or more than once on all three planes.
    float _all_planes_cosmic_pixels_tagged_once_or_more     = 0.;
    float _all_planes_cosmic_pixels_tagged_more_than_once   = 0.; 
    float _all_planes_neutrino_pixels_tagged_once_or_more   = 0.;
    float _all_planes_neutrino_pixels_tagged_more_than_once = 0.;
    float _all_planes_all_pixels_tagged_once_or_more        = 0.;
    float _all_planes_all_pixels_tagged_more_than_once      = 0.;

    // Start a loop over the pixels on each of the planes by looping over each of the planes.

    // Loop over the planes.
    for (size_t plane = 0; plane < 3; plane++) {

      // Declare three empty images: (1) Neutrino Pixel Image, (2) Cosmic Pixel Image (Based on which was tagged), (3) Entire Pixel Image.
      // I fill the image with the same dimensions as the image at this location in imgs_v.at(plane) ('plane' corresponds to the plane number).                   
      larcv::Image2D neutrino_pixels_tagged_info_img( imgs_v.at(plane).meta() );
      larcv::Image2D cosmic_pixels_tagged_info_img( imgs_v.at(plane).meta() );
      larcv::Image2D all_pixels_tagged_info_img( imgs_v.at(plane).meta() );

      // Empty each of these vectors.
      neutrino_pixels_tagged_info_img.paint(0.0);
      cosmic_pixels_tagged_info_img.paint(0.0);
      all_pixels_tagged_info_img.paint(0.0);

      // Declare values for the number of cosmic, neutrino, and all 'truth' pixels on the three planes.
      float _single_plane_cosmic_pixels          = 0.;
      float _single_plane_neutrino_pixels        = 0.;
      float _single_plane_all_pixels             = 0.;

      // Declare values for the number of cosmic, neutrino, and all 'tagged' pixels on the three planes.
      float _single_plane_cosmic_pixels_tagged   = 0.;
      float _single_plane_neutrino_pixels_tagged = 0.;
      float _single_plane_all_pixels_tagged      = 0.;

      // Declare variables for the number of pixels (of each type) tagged once or more than once on a single plane.
      float _single_plane_cosmic_pixels_tagged_once_or_more     = 0.;
      float _single_plane_cosmic_pixels_tagged_more_than_once   = 0.;
      float _single_plane_neutrino_pixels_tagged_once_or_more   = 0.;
      float _single_plane_neutrino_pixels_tagged_more_than_once = 0.;
      float _single_plane_all_pixels_tagged_once_or_more        = 0.;
      float _single_plane_all_pixels_tagged_more_than_once      = 0.;

      // Define five types of images on this plane: (1) the StopMu image, (2) the ThruMu image, (3) the DeadCh image, (4) the Seg image, and (5) the Event image.
      const larcv::Image2D& badchimg  = badch_v.at(plane);
      const larcv::Image2D& segimg    = segs_v.at(plane);
      const larcv::Image2D& eventimg  = imgs_v.at(plane);

      // First, deal with the truth info.  Count the number of cosmic pixels, neutrino pixels, and all pixels there are on the three planes.
      // This involves looping over the 'eventimg', seeing if a pixel is above threshold, seeing if it is a neutrino pixel or a cosmic pixel, and then incrementing the type of pixels for the correct category.
      // I do not have to worry about redundancy with the pixels here.
      
      // Loop over the rows.
      for (size_t event_row_iter = 0; event_row_iter < eventimg.meta().rows(); event_row_iter++) {

	// Loop over the columns.
	for (size_t event_col_iter = 0; event_col_iter < eventimg.meta().cols(); event_col_iter++) {

	  // Check to see if the pixel is above threshold, defined above from the config file as 'fthreshold', and if it is not a dead channel.
	  // If so, then continue on in the loop.
	  if (eventimg.pixel(event_row_iter, event_col_iter) < fthreshold || badchimg.pixel(event_row_iter, event_col_iter) > 0) { continue; }

	  // If this is the case, then the pixel is either a 'truth' cosmic pixel or a truth neutrino pixel.  I can increment the values of '_all_planes_all_pixels' and '_single_plane_all_pixels' for this reason.
	  _all_planes_all_pixels   += 1.0;
	  _single_plane_all_pixels += 1.0;

	  // Convert the value of 'event_row_iter' and 'event_col_iter' to 'x' and 'y' coordinates in the image.
	  float truth_pix_x = eventimg.meta().pos_x(event_col_iter);
	  float truth_pix_y = eventimg.meta().pos_y(event_row_iter);

	  // Check to see if it is a neutrino pixel using the 'isNeutrinoPixel' function.
	  bool is_truth_neutrino_pixel = isNeutrinoPixel(segimg, truth_pix_x, truth_pix_y);

	  // Using this information, you can increment the correct variables.
	  if (is_truth_neutrino_pixel == true) {

	    // Increment '_single_plane_neutrino_pixels' and '_all_planes_neutrino_pixels'.
	    _single_plane_neutrino_pixels += 1.0;
	    _all_planes_neutrino_pixels   += 1.0;

	  } // End of the conditional for if the pixel is a neutrino pixel.

	  if (is_truth_neutrino_pixel == false) {

	    // Increment '_single_plane_cosmic_pixels' and '_all_planes_cosmic_pixels'.
	    _single_plane_cosmic_pixels += 1.0;
	    _all_planes_cosmic_pixels   += 1.0;

	  } // End of the conditional for if the pixel is a cosmic pixel.

	} // End of the loop over the columns.

      } // End of the loop over the rows.

      // Finally declare the vector of 'Pixel2DCluster' objects, which themselves are vectors of type 'Pixel2D'.  These contain the row and column of the places in the image that the track passed through.

      // Use the object declared above to obtain the PixelCluster2D list from the 'EventPixel2D' class.  The name of that object is 'ev_pixel_clusters'.
      const std::vector<larcv::Pixel2DCluster>& ev_thrumu_pixel_info_from_event   =  ev_thrumu_pixel_clusters->Pixel2DClusterArray(plane);
      const std::vector<larcv::Pixel2DCluster>& ev_stopmu_pixel_info_from_event   =  ev_stopmu_pixel_clusters->Pixel2DClusterArray(plane);

      // Create a vector of these two clusters so that we can use a loop.
      const std::vector< std::vector < larcv::Pixel2DCluster > > plane_pixel_info_from_event_vector{ev_thrumu_pixel_info_from_event, ev_stopmu_pixel_info_from_event};

      for (size_t pixel_cluster_vector_counter = 0; pixel_cluster_vector_counter < 2; pixel_cluster_vector_counter++) {

	// Initialize 'plane_pixel_info_from_event' as the entry at index 'pixel_cluster_vector_counter' in  'plane_pixel_info_from_event_vector'.
	const std::vector<larcv::Pixel2DCluster>& plane_pixel_info_from_event = plane_pixel_info_from_event_vector.at(pixel_cluster_vector_counter);

	// Now, I will begin a loop over the tracks in order to fill the 'neutrino_pixels_tagged_info_img' and 'cosmic_pixels_tagged_info_img' vectors.
	for (size_t track_iter = 0; track_iter < plane_pixel_info_from_event.size(); track_iter++) {

	  // Declare an object of type 'Pixel2DCluster' for this entry of 'plane_pixel_info_from_event'.
	  const larcv::Pixel2DCluster& track_pixels = plane_pixel_info_from_event.at(track_iter);
	
	  // Start a loop over the pixels in this vector.
	  for (size_t pixel_iter = 0; pixel_iter < track_pixels.size(); pixel_iter++) {

	    // Declare a variable for the pixel at this particular value.
	    const larcv::Pixel2D& pixel = track_pixels.at(pixel_iter);

	    // Declare a variable for the row and the col for this particular pixel.
	    // col - x()
	    size_t col = pixel.X();
	    // row - y()
	    size_t row = pixel.Y();

	    // Convert these to the 'x' and 'y' coordinates in the image.
	    float x = eventimg.meta().pos_x(col);
	    float y = eventimg.meta().pos_y(row);

	    // Check if this pixel is a neutrino pixel. 
	    // Declare a boolean for if it is a neutrino pixel and find its value from the 'isNeutrinoPixel' function.
	    bool is_neutrino_pixel = isNeutrinoPixel(segimg, x, y);

	    // See if the pixel has been tagged if it is a neutrino pixel.
	    if (is_neutrino_pixel == true) {

	      // Check to make sure that the pixel is not a dead pixel.  Continue if it is.
	      if (fabs(badchimg.pixel(row,col)) > 0) { continue; }

	      // Increment '_all_planes_neutrino_pixels_tagged', '_single_plane_neutrino_pixels_tagged', '_all_planes_all_pixels_tagged', and '_single_plane_all_pixels_tagged' if this entry of the 'neutrino_pixels_tagged_info_img' is less than 1.0 and this is not a bad channel.
	      if (fabs(neutrino_pixels_tagged_info_img.pixel(row, col)) < 0.5) {

		_single_plane_neutrino_pixels_tagged += 1.0;
		_all_planes_neutrino_pixels_tagged   += 1.0;
		
		_single_plane_all_pixels_tagged      += 1.0;
		_all_planes_all_pixels_tagged        += 1.0;

	      } // End of the loop over the count of neutrino pixels.

	      // Find the value of the pixel at this entry:
	      _current_neutrino_pixel_value = neutrino_pixels_tagged_info_img.pixel(row, col);
	      _current_all_pixel_value      = all_pixels_tagged_info_img.pixel(row, col);

	      // Use the 'set_pixel' function to adjust the values of the 'neutrino_pixels_tagged_info_img' and 'all_pixels_tagged_info_img'.
	      neutrino_pixels_tagged_info_img.set_pixel(row, col, _current_neutrino_pixel_value + 1.0);
	      all_pixels_tagged_info_img.set_pixel(row, col, _current_all_pixel_value + 1.0);
	      
	    } // End of the conditional for if the pixel is a neutrino pixel.


	  // If 'is_neutrino_pixel' is false, then this is a cosmic pixel.  I can fill out the analogous information for the cosmic pixels.
	  // This is coded very literally so that I do not run into problems with equality operators.
	  if (is_neutrino_pixel == false) {

	    // Check to make sure that the pixel is not a dead pixel.  Continue if it is.
	    if (fabs(badchimg.pixel(row,col)) > 0) { continue; }

	    // Increment '_all_planes_cosmic_pixels_tagged', '_single_plane_cosmic_pixels_tagged', _all_planes_all_pixels_tagged', and '_single_plane_all_pixels_tagged' if this entry of the 'cosmic_pixels_tagged_info_img' is less than 1.0 and this is not a bad channel.
	    if (fabs(cosmic_pixels_tagged_info_img.pixel(row, col)) < 0.5) {

	      _single_plane_cosmic_pixels_tagged += 1.0;
	      _all_planes_cosmic_pixels_tagged   += 1.0;
	
	      _single_plane_all_pixels_tagged    += 1.0;
	      _all_planes_all_pixels_tagged      += 1.0;

	    } // End of the loop over the count of cosmic pixels.
	    
	    // Find the value of the pixel at this entry:                                                                                                                        
	    _current_cosmic_pixel_value   = cosmic_pixels_tagged_info_img.pixel(row, col);
	    _current_all_pixel_value      = all_pixels_tagged_info_img.pixel(row, col);
	    
	    // Use the 'set_pixel' function to adjust the values of the 'cosmic_pixels_tagged_info_img' and 'all_pixels_tagged_info_img'.      
	    cosmic_pixels_tagged_info_img.set_pixel(row, col, _current_cosmic_pixel_value + 1.0);
	    all_pixels_tagged_info_img.set_pixel(row, col, _current_all_pixel_value + 1.0);
	    
	  } // End of the conditional for if the pixel is a cosmic pixel.                                                                                                         
	  
	  } // End of the loop over the pixels in each of the Pixel2DCluster objects.
	  
	}  // End of the loops over the tracks within the event.

      } // End of the loop over both the thrumu and stopmu information for this event.

      // Now, I can do conditionals for which plane I am on and calculate the correct values based on that.

      // Use a loop over the 'cosmic_pixels_tagged_info_img' vector to count the number of pixels that are tagged once or more and those tagged more than once.
      for (size_t row_iter = 0; row_iter < cosmic_pixels_tagged_info_img.meta().rows(); row_iter++) {

	for (size_t col_iter = 0; col_iter < cosmic_pixels_tagged_info_img.meta().cols(); col_iter++) {

	  // Check to see if the pixel has been tagged at least once.
	  // '0.5' is the value that I will use to tell if a pixel has been tagged at least once.  It would be '0' if it has never been tagged.
	  if ( cosmic_pixels_tagged_info_img.pixel(row_iter, col_iter) > 0.5) {

	    // Increment '_single_plane_cosmic_pixels_tagged_once_or_more' and '_single_plane_all_pixels_tagged_once_or_more'.
	    _single_plane_cosmic_pixels_tagged_once_or_more += 1.0;
	    _single_plane_all_pixels_tagged_once_or_more    += 1.0;

	    // Also increment '_all_planes_cosmic_pixels_tagged_once_or_more' and '_all_planes_all_pixels_tagged_once_or_more'.  This will be a running variable.
	    _all_planes_cosmic_pixels_tagged_once_or_more   += 1.0;
	    _all_planes_all_pixels_tagged_once_or_more      += 1.0;

	  } // End of the conditional over cosmic pixels tagged once or more.

	  // Check to see if the pixel has been tagged more than once.
	  // '1.5' is the value that I will use to tell if a pixel has been tagged more than once.  It would be at least '2.0', but this value will give the same result.
	  if ( cosmic_pixels_tagged_info_img.pixel(row_iter, col_iter) > 1.5 ) {

	    // Increment '_single_plane_cosmic_pixels_tagged_more_than_once' and '_single_plane_all_pixels_tagged_more_than_once'.
	    _single_plane_cosmic_pixels_tagged_more_than_once += 1.0;
	    _single_plane_all_pixels_tagged_more_than_once    += 1.0;

	    // Also increment '_all_planes_cosmic_pixels_tagged_more_than_once' and '_all_planes_all_pixels_tagged_more_than_once'.  This will be a running variable.
	    _all_planes_cosmic_pixels_tagged_more_than_once += 1.0;
	    _all_planes_all_pixels_tagged_more_than_once    += 1.0;

	  } // End of the conditional over the cosmic pixels tagged more than once.

	} // End of the loop over the columns.

      } // End of the loop over the rows.
	
      // Use a loop over the 'neutrino_pixels_tagged_info_img' vector to count the number of pixels that are tagged once or more and those tagged more than once.
      // Use the same variables as before. (I have to redeclare them).
      for (size_t row_iter = 0; row_iter < neutrino_pixels_tagged_info_img.meta().rows(); row_iter++) {

	for (size_t col_iter = 0; col_iter < neutrino_pixels_tagged_info_img.meta().cols(); col_iter++) {

	  // Check to see if the pixel has been tagged more than once. 
	  // '0.5' is the value that I will use to tell if a pixel has been tagged at least once.  It would be '0' if it has never been tagged.
	  if ( neutrino_pixels_tagged_info_img.pixel(row_iter, col_iter) > 0.5 ) {

	    // Increment '_single_plane_neutrino_pixels_tagged_once_or_more' and '_single_plane_all_pixels_tagged_once_or_more'.                          
	    _single_plane_neutrino_pixels_tagged_once_or_more += 1.0;
	    _single_plane_all_pixels_tagged_once_or_more      += 1.0;

	    // Also increment '_all_planes_neutrino_pixels_tagged_once_or_more' and '_all_planes_all_pixels_tagged_once_or_more'.  This will be a running variable.    
	    _all_planes_neutrino_pixels_tagged_once_or_more   += 1.0;
	    _all_planes_all_pixels_tagged_once_or_more        += 1.0;
	      
	  } // End of the conditional over the neutrino pixels tagged once or more.

	  // Check to see if the pixel has been tagged more than once.                                                      
	  // '1.5' is the value that I will use to tell if a pixel has been tagged more than once.  It would be at least '2.0', but this value will give the same result.   
	  if ( neutrino_pixels_tagged_info_img.pixel(row_iter, col_iter) > 1.5 ) {

	    // Increment '_single_plane_neutrino_pixels_tagged_more_than_once'.                                         
	    _single_plane_neutrino_pixels_tagged_more_than_once += 1.0;
	    _single_plane_all_pixels_tagged_more_than_once      += 1.0;

	    // Also increment '_all_planes_neutrino_pixels_tagged_more_than_once' and '_all_planes_all_pixels_tagged_more_than_once'.  This will be a running variable.  
	    _all_planes_neutrino_pixels_tagged_more_than_once   += 1.0;
	    _all_planes_all_pixels_tagged_more_than_once        += 1.0;

	  } // End of the contional over the neutrino pixels that are tagged more than once.

	} // End of the loop over the columns.

      } // End of the loop over the rows.

	      
      // Calculate the correct fractions based on which plane I am on.  You can calculate all of the fractions from the variables that have been formed.
      if (plane == 0) {

	_u_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once/_single_plane_all_pixels;
	_u_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels;
	_u_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once/_single_plane_cosmic_pixels;
	_u_plane_frac_cosmic_pixels_not_tagged            = 1.0 - (_single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels);
	_u_plane_frac_neutrino_pixels_not_tagged          = 1.0 - (_single_plane_neutrino_pixels_tagged_once_or_more/_single_plane_neutrino_pixels);
	
	// Set the information for each type of pixel for the u-plane.
	_u_plane_all_pixels                               = _single_plane_all_pixels;
	_u_plane_cosmic_pixels                            = _single_plane_cosmic_pixels;
	_u_plane_neutrino_pixels                          = _single_plane_neutrino_pixels;
	 
      } // End of the conditional for the u-plane pixels.

      if (plane == 1) {

	_v_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once/_single_plane_all_pixels;
	_v_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels;
	_v_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once/_single_plane_cosmic_pixels;
	_v_plane_frac_cosmic_pixels_not_tagged            = 1.0 - (_single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels);
	_v_plane_frac_neutrino_pixels_not_tagged          = 1.0 - (_single_plane_neutrino_pixels_tagged_once_or_more/_single_plane_neutrino_pixels);

	// Set the information for each type of pixel for the v-plane.                                                               
	_v_plane_all_pixels                         = _single_plane_all_pixels;
	_v_plane_cosmic_pixels                      = _single_plane_cosmic_pixels;
	_v_plane_neutrino_pixels                    = _single_plane_neutrino_pixels;


      } // End of the conditional for the v-plane pixels.

      if (plane == 2) {

	_y_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once/_single_plane_all_pixels;
	_y_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels;
	_y_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once/_single_plane_cosmic_pixels;
	_y_plane_frac_cosmic_pixels_not_tagged            = 1.0 - (_single_plane_cosmic_pixels_tagged_once_or_more/_single_plane_cosmic_pixels);
	_y_plane_frac_neutrino_pixels_not_tagged          = 1.0 - (_single_plane_neutrino_pixels_tagged_once_or_more/_single_plane_neutrino_pixels);

	// Set the information for each type of pixel for the y-plane.                                                         
	_y_plane_all_pixels                         = _single_plane_all_pixels;
	_y_plane_cosmic_pixels                      = _single_plane_cosmic_pixels;
	_y_plane_neutrino_pixels                    = _single_plane_neutrino_pixels;


      } // End of the conditional for the y-plane pixels.

    } // End of the loop over the planes.

    // Calculate the same quantites for all of the planes using the '_all_planes' variables.
    _frac_all_pixels_tagged_more_than_once    = _all_planes_all_pixels_tagged_more_than_once/_all_planes_all_pixels;
    _frac_cosmic_pixels_tagged_once_or_more   = _all_planes_cosmic_pixels_tagged_once_or_more/_all_planes_cosmic_pixels;
    _frac_cosmic_pixels_tagged_more_than_once = _all_planes_cosmic_pixels_tagged_more_than_once/_all_planes_cosmic_pixels;
    _frac_cosmic_pixels_not_tagged            = 1.0 - (_all_planes_cosmic_pixels_tagged_once_or_more/_all_planes_cosmic_pixels);
    _frac_neutrino_pixels_not_tagged          = 1.0 - (_all_planes_neutrino_pixels_tagged_once_or_more/_all_planes_neutrino_pixels);

    _all_planes_all_pixels_tree               = _all_planes_all_pixels;
    _all_planes_cosmic_pixels_tree            = _all_planes_cosmic_pixels;
    _all_planes_neutrino_pixels_tree          = _all_planes_neutrino_pixels;

    // Print out the 'truth' number of the each type of pixels.
    std::cout << "Number of neutrino pixels = " << _all_planes_neutrino_pixels << "." << std::endl;
    std::cout << "Number of cosmic pixels = " << _all_planes_cosmic_pixels << "." << std::endl;
    std::cout << "Number of pixels altogether = " << _all_planes_all_pixels << "." << std::endl;
    std::cout << "\n" << std::endl;

    // Print out the 'tagged' number of each type of pixels.
    std::cout << "Number of tagged neutrino pixels = " << _all_planes_neutrino_pixels_tagged << "." << std::endl;
    std::cout << "Number of tagged cosmic pixels = " << _all_planes_cosmic_pixels_tagged <<"." << std::endl;
    std::cout << "Number of tagged pixels altogether = " << _all_planes_all_pixels_tagged << "." << std::endl;
    std::cout << "\n" << std::endl;

    tree->Fill();

    if ( ientry>=100 )
      break;
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main

// Define a function which will test if a pixel is a neutrino pixel, called 'isNeutrinoPixel'.                                
// Input parameters: 'segmented_image' - type 'larcv::Image2D', this is the image that contains only the neutrino pixels.                                  
//                   'pix_x'           - type 'float', this is the x-coordinate of the pixel of which you are determining the type.          
//                   'pix_y'           - type 'float', this is the y-coordinate of the pixel of which you are determining the type.                                      
bool isNeutrinoPixel(const larcv::Image2D& segmented_image, float pix_x, float pix_y) {

  // Declare a boolean for if this a neutrino pixel and initialize it to 'false'.                          
  bool is_neutrino_pixel = false;

  // Declare a boolean for if the pixel is in the 'seg_image'.                                                
  bool in_seg_image = false;
  int seg_row = -1;
  int seg_col = -1;

  if ( pix_x>segmented_image.meta().min_x() && pix_x<segmented_image.meta().max_x()
       && pix_y>segmented_image.meta().min_y() && pix_y<segmented_image.meta().max_y() ) {
    in_seg_image = true;

    // Find the coordinates in the 'seg_image' now.                                                                                
    seg_row = segmented_image.meta().row(pix_y);
    seg_col = segmented_image.meta().col(pix_x);
  }

  // If the particle is a neutrino pixel, then you can set 'is_neutrino_pixel' to 'true'.                         
  if ( in_seg_image == true  && segmented_image.pixel(seg_row,seg_col) > 0 ) {

    is_neutrino_pixel = true;

  }

  return is_neutrino_pixel;

}


