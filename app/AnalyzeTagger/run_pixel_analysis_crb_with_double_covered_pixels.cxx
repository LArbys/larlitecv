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

  // run pixel analysis. use 

  // we need to have a data coordinator for each stage because the number of entries could be different.
  // we'll coordinate by using event,subrun,run information
  std::string data_folder = "~/working/data/larbys/cosmic_tagger_dev/";    // blade
  //std::string data_folder = "~/data/larbys/cosmic_tagger/mcc7_bnbcosmic/"; // nudot

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
  dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv_p1.root", "larcv" ); //stopmu-tagger output

  larlitecv::DataCoordinator dataco_nucand;
  dataco_nucand.add_inputfile( "../ContainedROI/bin/output_containedroi_larcv.root", "larcv");
  dataco_nucand.add_inputfile( "../ContainedROI/bin/output_containedroi_larlite.root", "larlite");

  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_thrumu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_stopmu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_nucand.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();
  dataco_nucand.initialize();

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;
  std::cout << "data[nucand] entries=" << dataco_nucand.get_nentries("larcv") << std::endl;

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet pixana_cfg = cfg.get<larcv::PSet>("PixelAnalysis");
  float fthreshold   = pixana_cfg.get<float>("PixelThreshold");
  int fvertex_radius = pixana_cfg.get<int>("PixelRadius");
  int verbosity      = pixana_cfg.get<int>("Verbosity",0);

  // setup output

  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
  TTree* tree = new TTree("pixana", "Pixel-level analysis");
  int run, subrun, event;
  int lepton_boundary; // index of boundary that lepton end point is nearest
  int vertex_boundary; // index of boundary that vertex is nearest
  int ncosmic_pixels;  // number of non-neutrino pixels
  int nnu_pixels;      // number of neutrino pixels
  int nvertex_pixels;  // number of neutrino pixels within some pixel radius of vertex
  int ncosmic_tagged;  // number of non-neutrino pixels tagged
  int nnu_tagged;      // number of neutrino pixels tagged
  int nvertex_tagged;  // number of neutrino pixels within some pixel radius of vertex tagged
  int mode;            // interaction mode
  int current;         // interaction cufrrent
  int num_protons_over60mev; // as named
  int num_rois;        // number of identified ROis
  float EnuGeV;        // neutrino energy in GeV
  float fdwall;        // dwall
  float dwall_lepton;  // dwll for end of lepton
  float frac_cosmic;   // fraction of cosmic tagged
  float frac_nu;       // fraction of neutrino pixels tagged
  float frac_vertex;   // fraction of vertex pixels tagged
  float frac_inroi;    // fraction of pixels contained in one of the rois
  float primary_proton_ke; // ke of the proton driving from the hit nucleon
  float lepton_cosz;
  float lepton_phiz;
  float fpos[3];       // vertex

  // Declare variables for the quantities that I am interested in plotting.  These are for the quantities including ALL of the planes.
  float _frac_all_pixels_tagged_more_than_once;
  float _frac_cosmic_pixels_tagged_once_or_more;
  float _frac_cosmic_pixels_tagged_more_than_once;
  float _frac_cosmic_pixels_not_tagged;
  float _frac_neutrino_pixels_not_tagged;

  // Declare the same information but specific to each plane.
  // u
  float _u_plane_frac_all_pixels_tagged_more_than_once;
  float _u_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _u_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _u_plane_frac_cosmic_pixels_not_tagged;
  float _u_plane_frac_neutrino_pixels_not_tagged;

  // v
  float _v_plane_frac_all_pixels_tagged_more_than_once;
  float _v_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _v_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _v_plane_frac_cosmic_pixels_not_tagged;
  float _v_plane_frac_neutrino_pixels_not_tagged;

  // y
  float _y_plane_frac_all_pixels_tagged_more_than_once;
  float _y_plane_frac_cosmic_pixels_tagged_once_or_more;
  float _y_plane_frac_cosmic_pixels_tagged_more_than_once;
  float _y_plane_frac_cosmic_pixels_not_tagged;
  float _y_plane_frac_neutrino_pixels_not_tagged;

  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");
  tree->Branch("ncosmic_pixels",&ncosmic_pixels,"ncosmic_pixels/I");
  tree->Branch("nnu_pixels",&nnu_pixels,"nnu_pixels/I");
  tree->Branch("nvertex_pixels",&nvertex_pixels,"nvertex_pixels/I");
  tree->Branch("ncosmic_tagged",&ncosmic_tagged,"ncosmic_tagged/I");
  tree->Branch("nnu_tagged",&nnu_tagged,"nnu_tagged/I");
  tree->Branch("nvertex_tagged",&nvertex_tagged,"nvertex_tagged/I");
  tree->Branch("num_rois", &num_rois, "num_rois/I");
  tree->Branch("mode",&mode,"mode/I");
  tree->Branch("current",&current,"current/I");
  tree->Branch("lepton_boundary",&lepton_boundary,"lepton_boundary/I");
  tree->Branch("vertex_boundary",&vertex_boundary,"vertex_boundary/I");
  tree->Branch("num_protons_over60mev", &num_protons_over60mev, "num_protons_over60mev/I");
  tree->Branch("EnuGeV",&EnuGeV,"EnuGeV/F");
  tree->Branch("dwall",&fdwall,"dwall/F");
  tree->Branch("dwall_lepton",&dwall_lepton,"dwall_lepton/F");
  tree->Branch("frac_cosmic",&frac_cosmic,"frac_cosmic/F");
  tree->Branch("frac_nu",&frac_nu,"frac_nu/F");
  tree->Branch("frac_vertex",&frac_vertex,"frac_vertex/F");
  tree->Branch("frac_inroi", &frac_inroi, "frac_inroi/F");
  tree->Branch("primary_proton_ke", &primary_proton_ke, "primary_proton_ke/F");
  tree->Branch("pos",fpos,"pos[3]/F");
  tree->Branch("lepton_cosz", &lepton_cosz, "lepton_cosz/F");
  tree->Branch("lepton_phiz", &lepton_phiz, "lepton_phiz/F");

  // Create branches for the quantities that I am interested in plotting.
  tree->Branch("_frac_all_pixels_tagged_more_than_once", &_frac_all_pixels_tagged_more_than_once, "frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_frac_cosmic_pixels_tagged_once_or_more", &_frac_cosmic_pixels_tagged_once_or_more, "frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_frac_cosmic_pixels_tagged_more_than_once", &_frac_cosmic_pixels_tagged_more_than_once, "frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_frac_cosmic_pixels_not_tagged", &_frac_cosmic_pixels_not_tagged, "frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_frac_neutrino_pixels_not_tagged", &_frac_neutrino_pixels_not_tagged, "frac_neutrino_pixels_not_tagged/F");

  // Create branches for the same quantities but for each of the wire planes.
  // u
  tree->Branch("_u_plane_frac_all_pixels_tagged_more_than_once", &_u_plane_frac_all_pixels_tagged_more_than_once, "u_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_tagged_once_or_more", &_u_plane_frac_cosmic_pixels_tagged_once_or_more, "u_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_tagged_more_than_once", &_u_plane_frac_cosmic_pixels_tagged_more_than_once, "u_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_u_plane_frac_cosmic_pixels_not_tagged", &_u_plane_frac_cosmic_pixels_not_tagged, "u_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_u_plane_frac_neutrino_pixels_not_tagged", &_u_plane_frac_neutrino_pixels_not_tagged, "u_plane_frac_neutrino_pixels_not_tagged/F");

  // v
  tree->Branch("_v_plane_frac_all_pixels_tagged_more_than_once", &_v_plane_frac_all_pixels_tagged_more_than_once, "v_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_tagged_once_or_more", &_v_plane_frac_cosmic_pixels_tagged_once_or_more, "v_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_tagged_more_than_once", &_v_plane_frac_cosmic_pixels_tagged_more_than_once, "v_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_v_plane_frac_cosmic_pixels_not_tagged", &_v_plane_frac_cosmic_pixels_not_tagged, "v_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_v_plane_frac_neutrino_pixels_not_tagged", &_v_plane_frac_neutrino_pixels_not_tagged, "v_plane_frac_neutrino_pixels_not_tagged/F");

  // y
  tree->Branch("_y_plane_frac_all_pixels_tagged_more_than_once", &_y_plane_frac_all_pixels_tagged_more_than_once, "y_plane_frac_all_pixels_tagged_more_than_once/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_tagged_once_or_more", &_y_plane_frac_cosmic_pixels_tagged_once_or_more, "y_plane_frac_cosmic_pixels_tagged_once_or_more/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_tagged_more_than_once", &_y_plane_frac_cosmic_pixels_tagged_more_than_once, "y_plane_frac_cosmic_pixels_tagged_more_than_once/F");
  tree->Branch("_y_plane_frac_cosmic_pixels_not_tagged", &_y_plane_frac_cosmic_pixels_not_tagged, "y_plane_frac_cosmic_pixels_not_tagged/F");
  tree->Branch("_y_plane_frac_neutrino_pixels_not_tagged", &_y_plane_frac_neutrino_pixels_not_tagged, "y_plane_frac_neutrino_pixels_not_tagged/F");


  int nentries = dataco_stopmu.get_nentries("larcv");

  dataco_source.goto_entry(0,"larcv");
  dataco_thrumu.goto_entry(0,"larcv");
  dataco_stopmu.goto_entry(0,"larcv");
  dataco_nucand.goto_entry(0,"larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco_nucand.goto_entry(ientry,"larcv");

    dataco_nucand.get_id(run,subrun,event);

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << std::endl;
      std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    dataco_thrumu.goto_event(run,subrun,event,"larcv");
    dataco_source.goto_event(run,subrun,event,"larcv");
    dataco_nucand.goto_event(run,subrun,event,"larcv");

    // initialize the output variables
    ncosmic_pixels = 0;
    nnu_pixels = 0;
    nvertex_pixels = 0;
    ncosmic_tagged = 0;
    nnu_tagged = 0;
    nvertex_tagged = 0;
    mode = 0;
    current = 0;
    EnuGeV = 0.;
    fdwall = 0.;
    frac_cosmic = 0.;
    frac_vertex = 0.;
    frac_nu = 0.;
    fpos[0] = fpos[1] = fpos[2] = 0.;
    num_protons_over60mev = 0;
    primary_proton_ke = 0.;

    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_segs   = (larcv::EventImage2D*)dataco_source.get_larcv_data(larcv::kProductImage2D,"seg_comb_tpc"); // Previously 'segment' in the second entry.
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"modimgs");
    larcv::EventImage2D* ev_badch  = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"badch");	
    larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"marked3d");
    larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D,"stopmu");  	

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>& segs_v   = ev_segs->Image2DArray();
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
    const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();

    // Declare variables for the total number of both types of pixels, and pixels in general, on each of the planes.
    // These data types will have to be of type 'float' so that they can be used in the fractional values that I am looking at.
    float _all_planes_cosmic_pixels   = 0.;
    float _all_planes_neutrino_pixels = 0.;
    float _all_planes_all_pixels      = 0.;

    // Declare variables for the number of pixels (of each type) tagged once or more than once on all three planes.
    float _all_planes_cosmic_pixels_tagged_once_or_more     = 0.;
    float _all_planes_cosmic_pixels_tagged_more_than_once   = 0.; 
    float _all_planes_neutrino_pixels_tagged_once_or_more   = 0.;
    float _all_planes_neutrino_pixels_tagged_more_than_once = 0.;
    float _all_planes_all_pixels_tagged_once_or_more        = 0.;
    float _all_planes_all_pixels_tagged_more_than_once      = 0.;

    // Declare variables for the number of pixels not_tagged on all of the planes.
    float _all_planes_cosmic_pixels_not_tagged   = 0.;
    float _all_planes_neutrino_pixels_not_tagged = 0.;
    float _all_planes_all_pixels_not_tagged      = 0.;

    // Start a loop over the pixels on each of the planes by looping over each of the planes.

    // Loop over the planes.
    for (size_t plane = 0; plane < 2; plane++) {

      // Declare three empty images: (1) Neutrino Pixel Image, (2) Cosmic Pixel Image (Based on which was tagged), (3) Entire Pixel Image.
      // I fill the image with the same dimensions as the image at this location in imgs_v.at(plane) ('plane' corresponds to the plane number).                   
      larcv::Image2D neutrino_pixels_tagged_info_img( imgs_v.at(plane).meta() );
      larcv::Image2D cosmic_pixels_tagged_info_img( imgs_v.at(plane).meta() );
      larcv::Image2D all_pixels_tagged_info_img( imgs_v.at(plane).meta() );

      // Empty each of these vectors.
      neutrino_pixels_tagged_info_img.paint(0.0);
      cosmic_pixels_tagged_info_img.paint(0.0);
      all_pixels_tagged_info_img.paint(0.0);

      // Declare values for the number of cosmic, neutrino, and all pixels on the three planes.
      float _single_plane_cosmic_pixels   = 0.;
      float _single_plane_neutrino_pixels = 0.;
      float _single_plane_all_pixels      = 0.;

      // Declare variables for the number of pixels (of each type) tagged once or more than once on a single plane.
      float _single_plane_cosmic_pixels_tagged_once_or_more     = 0.;
      float _single_plane_cosmic_pixels_tagged_more_than_once   = 0.;
      float _single_plane_neutrino_pixels_tagged_once_or_more   = 0.;
      float _single_plane_neutrino_pixels_tagged_more_than_once = 0.;
      float _single_plane_all_pixels_tagged_once_or_more        = 0.;
      float _single_plane_all_pixels_tagged_more_than_once      = 0.;

      // Declare variables for the number of pixels not_tagged on all of the planes.                                                                
      float _single_plane_cosmic_pixels_not_tagged   = 0.;
      float _single_plane_neutrino_pixels_not_tagged = 0.;
      float _single_plane_all_pixels_not_tagged      = 0.;

      // Define five types of images on this plane: (1) the StopMu image, (2) the ThruMu image, (3) the DeadCh image, (4) the Seg image, and (5) the Event image.
      const larcv::Image2D& stopmuimg = stopmu_v.at(plane);
      const larcv::Image2D& thrumuimg = thrumu_v.at(plane);
      const larcv::Image2D& badchimg  = badch_v.at(plane);
      const larcv::Image2D& segimg    = segs_v.at(plane);
      const larcv::Image2D& eventimg  = imgs_v.at(plane);

      // Set the plane value to the correct type based on which plane we're on.
      const ::larcv::PlaneID_t plane_correct_type;

      if (plane == 0) {

	plane_correct_type = 0;

      }

      if (plane == 1) {

	plane_correct_type = 1;

      } 

      if (plane == 2) {

	plane_correct_type = 2;

      } 

      // Finally declare the vector of 'Pixel2DCluster' objects, which themselves are vectors of type 'Pixel2D'.  These contain the row and column of the places in the image that the track passed through.
      const std::vector<larcv::Pixel2DCluster>& plane_pixel_info_from_event   = const std::vector<larcv::Pixel2DCluster>& Pixel2DClusterArray(plane_correct_type);

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
	  size_t col = pixel.y();
	  // row - y()
	  size_t row = pixel.x();

	  // Check if this pixel is a neutrino pixel. 
	  // Declare a boolean for if it is a neutrino pixel.
	  bool is_neutrino_pixel = false;

	  // Find the x-coordinate and y-coordinate in the event image, and then convert it into coordinates in the 'segimg'.
	  float x = eventimg.meta().pos_x(col);
          float y = eventimg.meta().pos_y(row);
          bool in_seg_image = false;
          int seg_row = -1;
          int seg_col = -1;
          if ( x>segimg.meta().min_x() && x<segimg.meta().max_x()
	       && y>segimg.meta().min_y() && y<segimg.meta().max_y() ) {
            in_seg_image = true;

	    // Find the coordinates in the 'segimg' now.
            seg_row = segimg.meta().row(y);
            seg_col = segimg.meta().col(x);
          }
	  
	  // If the particle is a neutrino pixel, then you can set 'is_neutrino_pixel' to 'true'. 
          if ( in_seg_image && segimg.pixel(seg_row,seg_col)>0 ) {

	    is_neutrino_pixel = true;

	  }

	  // See if the pixel has been tagged if it is a neutrino pixel.
	  if (is_neutrino_pixel) {

	    // Check to make sure that the pixel is not a dead pixel.
	    if (fabs(badchimg.pixel(row,col)) < 0.001) {

	      // Increment '_all_planes_neutrino_pixels', '_single_plane_neutrino_pixels', '_all_planes_all_pixels', and '_single_plane_all_pixels' if this entry of the 'neutrino_pixels_tagged_info_img' is approximately 0 and this is not a bad channel.
	      if (fabs(neutrino_pixels_tagged_info_img.pixel(row, col)) < 0.001) {

		_single_plane_neutrino_pixels += 1.0;
		_all_planes_neutrino_pixels   += 1.0;
		
		_single_plane_all_pixels      += 1.0;
		_all_planes_all_pixels        += 1.0;

		// Check the other case now in which they both are less than 0.001, which means that this is a not_tagged neutrino pixel.                             
		if (stopmuimg.pixel(row,col) < 0.001 && thrumuimg.pixel(row,col) < 0.001) {

		  _single_plane_neutrino_pixels_not_tagged += 1.0;
		  _all_planes_neutrino_pixels_not_tagged   += 1.0;

		  _single_plane_all_pixels_not_tagged      += 1.0;
		  _all_planes_all_pixels_not_tagged        += 1.0;

		} // End of the loop over the untagged pixels.                                                                                                                      

	      } // End of the loop over the count of neutrino pixels.

	      // The values before had to be greater than 0 to see if the pixel had been tagged.  But for now, they just have to be greater than 0.001.
	      if (stopmuimg.pixel(row,col) > 0.001  || thrumuimg.pixel(row,col) > 0.001) {

		neutrino_pixels_tagged_info_img.pixel(row, col) += 1;
		all_pixels_tagged_info_img.pixel(row, col)      += 1;

	      } // End of the loop over if the pixel is tagged by ThruMu or StopMu.

	    } // End of the conditional for if the pixel is dead.

	  } // End of the conditional for if the pixel is a neutrino pixel.


	  // If 'is_neutrino_pixel' is false, then this is a cosmic pixel.  I can fill out the analogous information for the cosmic pixels.
	  if (!is_neutrino_pixel) {

	    // Check to make sure that the pixel is not a dead pixel.
	    if (fabs(badchimg.pixel(row,col)) < 0.001) {

	      // Increment '_all_planes_cosmic_pixels', '_single_plane_cosmic_pixels', _all_planes_all_pixels', and '_single_plane_all_pixels' if this entry of the 'cosmic_pixels_tagged_info_img' is approximately 0 and this is not a bad channel.
	      if (fabs(cosmic_pixels_tagged_info_img.pixel(row, col)) < 0.001) {

		_single_plane_cosmic_pixels += 1.0;
		_all_planes_cosmic_pixels   += 1.0;
	
		_single_plane_all_pixels    += 1.0;
		_all_planes_all_pixels      += 1.0;

		// Check the other case now in which they both are less than 0.001, which means that this is a not_tagged neutrino pixel.                                
		if (stopmuimg.pixel(row,col) < 0.001 && thrumuimg.pixel(row,col) < 0.001) {

                  _single_plane_cosmic_pixels_not_tagged += 1.0;
		  _all_planes_cosmic_pixels_not_tagged   += 1.0;

		  _single_plane_all_pixels_not_tagged    += 1.0;
		  _all_planes_all_pixels_not_tagged      += 1.0;

                } // End of the loop over the untagged pixels.                                                                                                                      

	      } // End of the loop over the count of cosmic pixels.
	     
	      // The values before had to be greater than 0 to see if the pixel had been tagged.  But for now, they just have to be greater than 0.001.                
	      if (stopmuimg.pixel(row,col) > 0.001  || thrumuimg.pixel(row,col) > 0.001) {

                cosmic_pixels_tagged_info_img.pixel(row, col)   += 1;
                all_pixels_tagged_info_img.pixel(row, col)      += 1;

              } // End of the loop over if the pixel is tagged by ThruMu or StopMu.                                                                                          
	    
	    } // End of the conditional for if the pixel is dead.                                                                                                            
	  
	  } // End of the conditional for if the pixel is a cosmic pixel.                                                                                                         
	  
	} // End of the loop over the pixels in each of the Pixel2DCluster objects.

      } // End of the loops over the tracks within the event.

      // Now, I can do conditionals for which plane I am on and calculate the correct values based on that.

	// Use a loop over the 'cosmic_pixels_tagged_info_img' vector to count the number of pixels that are tagged once or more and those tagged more than once.
	for (size_t row_iter = 0; row_iter < cosmic_pixels_tagged_info_img.meta().rows(); row_iter++) {

	  for (size_t col_iter = 0; col_iter < cosmic_pixels_tagged_info_img.meta().cols(); col_iter++) {

	    // Check to see if the pixel has been tagged at least once.
	    // '0.5' is the value that I will use to tell if a pixel has been tagged at least once.  It would be '0' if it has never been tagged.
	    if ( cosmic_pixels_tagged_info_img.pixel(row, col) > 0.5) {

	      // Increment '_single_plane_cosmic_pixels_tagged_once_or_more' and '_single_plane_all_pixels_tagged_once_or_more'.
	      _single_plane_cosmic_pixels_tagged_once_or_more += 1.0;
	      _single_plane_all_pixels_tagged_once_or_more    += 1.0;

	      // Also increment '_all_planes_cosmic_pixels_tagged_once_or_more' and '_all_planes_all_pixels_tagged_once_or_more'.  This will be a running variable.
	      _all_planes_cosmic_pixels_tagged_once_or_more   += 1.0;
	      _all_planes_all_pixels_tagged_once_or_more      += 1.0;

	    } // End of the conditional over cosmic pixels tagged once or more.

	    // Check to see if the pixel has been tagged more than once.
	    // '1.5' is the value that I will use to tell if a pixel has been tagged more than once.  It would be at least '2.0', but this value will give the same result.
	    else if ( cosmic_pixels_tagged_info_img.pixel(row, col) > 1.5 ) {

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
	// Use the same variables as before.
	for (row_iter = 0; row_iter < neutrino_pixels_tagged_info_img.meta().rows(); row_iter++) {

	  for (col_iter = 0; col_iter < neutrino_pixels_tagged_info_img.meta().cols(); col_iter++) {

	    // Check to see if the pixel has been tagged more than once. 
	    // '0.5' is the value that I will use to tell if a pixel has been tagged at least once.  It would be '0' if it has never been tagged.
	    if ( neutrino_pixels_tagged_info_img.pixel(row, col) > 0.5 ) {

	      // Increment '_single_plane_neutrino_pixels_tagged_once_or_more' and '_single_plane_all_pixels_tagged_once_or_more'.                          
	      _single_plane_neutrino_pixels_tagged_once_or_more += 1.0;
	      _single_plane_all_pixels_tagged_once_or_more      += 1.0;

	      // Also increment '_all_planes_neutrino_pixels_tagged_once_or_more' and '_all_planes_all_pixels_tagged_once_or_more'.  This will be a running variable.    
	      _all_planes_neutrino_pixels_tagged_once_or_more += 1.0;
	      _all_planes_all_pixels_tagged_once_or_more      += 1.0;
	      
	    } // End of the conditional over the neutrino pixels tagged once or more.

	    // Check to see if the pixel has been tagged more than once.                                                      
	    // '1.5' is the value that I will use to tell if a pixel has been tagged more than once.  It would be at least '2.0', but this value will give the same result.   
	    else if ( neutrino_pixels_tagged_info_img.pixel(row, col) > 1.5 ) {

	      // Increment '_single_plane_neutrino_pixels_tagged_more_than_once'                                         
	      _single_plane_neutrino_pixels_tagged_more_than_once += 1.0;
	      _single_plane_all_pixels_tagged_more_than_once    += 1.0;

	      // Also increment '_all_planes_neutrino_pixels_tagged_more_than_once' and '_all_planes_all_pixels_tagged_more_than_once'.  This will be a running variable.  
	      _all_planes_neutrino_pixels_tagged_more_than_once += 1.0;
	      _all_planes_all_pixels_tagged_more_than_once      += 1.0;

	      
	    } // End of the contional over the neutrino pixels that are tagged more than once.

	  } // End of the loop over the columns.

	} // End of the loop over the rows.

	      
	// Calculate the correct fractions based on which plane I am on.  You can calculate all of the fractions from the variables that have been formed.
	if (plane == 0) {

	  _u_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once    / _single_plane_all_pixels;
	  _u_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more   / _single_plane_cosmic_pixels;
	  _u_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once / _single_plane_cosmic_pixels;
	  _u_plane_frac_cosmic_pixels_not_tagged            = _single_plane_cosmic_pixels_not_tagged            / _single_plane_cosmic_pixels;
	  _u_plane_frac_neutrino_pixels_not_tagged          = _single_plane_neutrino_pixels_not_tagged          / _single_plane_neutrino_pixels;
	 
	} // End of the conditional for the u-plane pixels.

	if (plane == 1) {

          _v_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once    / _single_plane_all_pixels;
	  _v_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more   / _single_plane_cosmic_pixels;
          _v_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once / _single_plane_cosmic_pixels;
          _v_plane_frac_cosmic_pixels_not_tagged            = _single_plane_cosmic_pixels_not_tagged            / _single_plane_cosmic_pixels;
          _v_plane_frac_neutrino_pixels_not_tagged          = _single_plane_neutrino_pixels_not_tagged          / _single_plane_neutrino_pixels;

	} // End of the conditional for the v-plane pixels.

	if (plane == 2) {

          _y_plane_frac_all_pixels_tagged_more_than_once    = _single_plane_all_pixels_tagged_more_than_once    / _single_plane_all_pixels;
	  _y_plane_frac_cosmic_pixels_tagged_once_or_more   = _single_plane_cosmic_pixels_tagged_once_or_more   / _single_plane_cosmic_pixels;
          _y_plane_frac_cosmic_pixels_tagged_more_than_once = _single_plane_cosmic_pixels_tagged_more_than_once / _single_plane_cosmic_pixels;
          _y_plane_frac_cosmic_pixels_not_tagged            = _single_plane_cosmic_pixels_not_tagged            / _single_plane_cosmic_pixels;
          _y_plane_frac_neutrino_pixels_not_tagged          = _single_plane_neutrino_pixels_not_tagged          / _single_plane_neutrino_pixels;

	} // End of the conditional for the y-plane pixels.

    }

    // Calculate the same quantites for all of the planes using the '_all_planes' variables.
    _frac_all_pixels_tagged_more_than_once    = _all_planes_all_pixels_tagged_more_than_once    / _all_planes_all_pixels;
    _frac_cosmic_pixels_tagged_once_or_more   = _all_planes_cosmic_pixels_tagged_once_or_more   / _all_planes_cosmic_pixels;
    _frac_cosmic_pixels_tagged_more_than_once = _all_planes_cosmic_pixels_tagged_more_than_once / _all_planes_cosmic_pixels;
    _frac_cosmic_pixels_not_tagged            = _all_planes_cosmic_pixels_not_tagged            / _all_planes_cosmic_pixels;
    _frac_neutrino_pixels_not_tagged          = _all_planes_neutrino_pixels_not_tagged          / _all_planes_neutrino_pixels;


    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = (larcv::EventROI*)dataco_nucand.get_larcv_data(larcv::kProductROI,"containedroi");
    const std::vector<larcv::ROI>& containedrois_v = ev_contained_roi->ROIArray();
    num_rois = (int)containedrois_v.size();
    std::cout << "==ROIs==" << std::endl;
    for ( auto& roi : containedrois_v ) {
      std::cout << " roi: " << roi.dump();
    }

    // get other information, e.g. truth
    larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)dataco_source.get_larlite_data(larlite::data::kMCTruth,"generator");
    larlite::event_mctrack* ev_mctrack = (larlite::event_mctrack*)dataco_source.get_larlite_data(larlite::data::kMCTrack,"mcreco");

    // extract the truth quantities of interest
    const larlite::mcnu& neutrino = ev_mctruth->at(0).GetNeutrino();

    mode = neutrino.InteractionType();
    current = neutrino.CCNC();
    EnuGeV = neutrino.Nu().Momentum(0).E();
    const TLorentzVector& nu_pos = neutrino.Nu().Position();
    std::vector<float> fpos_v(3);
    std::vector<double> dpos(3);
    fpos_v[0] = nu_pos.X();
    fpos_v[1] = nu_pos.Y();
    fpos_v[2] = nu_pos.Z();
    dpos[0] = nu_pos.X();
    dpos[1] = nu_pos.Y();
    dpos[2] = nu_pos.Z();
    fdwall = dwall(fpos_v, vertex_boundary);
    if ( verbosity>0 )
      std::cout << " Enu=" << EnuGeV << std::endl;

    // get the vertex in the pixel coordinates
    std::vector<int> wid(3,-1);
    std::vector<int> vertex_col(3,-1);
    for (size_t p=0; p<3; p++) {
      wid[p] = ::larutil::Geometry::GetME()->WireCoordinate( dpos, p );
      if ( wid[p]>=0 && wid[p]<3456 )
	vertex_col[p] = imgs_v.at(p).meta().col(wid[p]);
      fpos[p] = fpos_v[p];
    }
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    float ftick = nu_pos[0]/cm_per_tick + 3200.0;
    int vertex_row = -1;
    if ( ftick >= imgs_v.at(0).meta().min_y() && ftick<=imgs_v.at(0).meta().max_y() )
      vertex_row = imgs_v.at(0).meta().row( ftick );

    // get the initial direction of the leptop
    const std::vector<larlite::mcpart>& particles = ev_mctruth->at(0).GetParticles();
    bool found_lepton = false;
    int hit_nucleon_id = 2;
    int lepton_track_id = -1;
    std::set<int> protonids;
    for ( auto  const& particle : particles ) {
      float KE = particle.Trajectory().front().E() - particle.Mass();
      if ( !found_lepton && (particle.PdgCode()==13 || particle.PdgCode()==-13) ) {
	// found the lepton
        const larlite::mctrajectory& traj = particle.Trajectory();
	std::cout << "  lepton E=" << particle.Trajectory().front().E() << " KE=" << KE << std::endl;
        found_lepton = true;
        lepton_track_id = particle.TrackId();
      }
      else if ( particle.PdgCode()==2212 ) {
        std::cout << "  a proton. p=" << particle.Momentum(0).Vect().Mag() 
          << " E=" << particle.Trajectory().front().E() << " KE=" << KE
          << " status=" << particle.StatusCode() 
          << " trackid=" << particle.TrackId() << " mother=" << particle.Mother()
          << std::endl;
        if ( particle.StatusCode()!=11 && KE>0.060 && protonids.find(particle.Mother())==protonids.end() ) 
	  num_protons_over60mev++; // status 11 means from genie? threshold cut. check that it isn't from the original proton
        protonids.insert( particle.TrackId() );
      }
      else if ( particle.PdgCode()==14 || particle.PdgCode()==-14 ) {
        std::cout << "  the neutrino (pdg=" << particle.PdgCode() << ") Enu=" << particle.Trajectory().front().E() << std::endl;
      } 
      else {
        std::cout << "  pdg=" << particle.PdgCode() 
          << " E=" << particle.Trajectory().front().E() << " KE=" << KE
          << " status=" << particle.StatusCode() 
          << " end process=" << particle.EndProcess()
          << " trackid=" << particle.TrackId() << " mother=" << particle.Mother() << std::endl;
      }

      // stuff we are saving
      if ( (particle.PdgCode()==2212 || particle.PdgCode()==2112) && particle.StatusCode()==11 ) {
        hit_nucleon_id = particle.TrackId();
      }
      if ( particle.PdgCode()==2212 && particle.Mother()==hit_nucleon_id ) {
        primary_proton_ke = KE;
        protonids.insert(particle.TrackId());
      }

    }//end of particle track loop

    std::cout << "lepton track id = " << lepton_track_id << std::endl;
    std::cout << "num_protons_over60mev=" << num_protons_over60mev << std::endl;
    std::cout << "primary_proton_ke=" << primary_proton_ke << std::endl;
    
    // loop over MC tracks, find the neutrino lepton by matching vertex
    for ( auto const& track : *ev_mctrack ) {
      if ( std::abs(track.PdgCode())!=13  ) continue;
      if ( track.size()==0 ) continue;
      const TLorentzVector& track_start = track.front().Position();
      std::vector<float> fstart(3);
      fstart[0] = track_start.X();
      fstart[1] = track_start.Y();
      fstart[2] = track_start.Z();

      float vert_dist = 0.;
      for (int v=0; v<3; v++) {
        float dv = fpos_v[v]-fstart[v];
        vert_dist += dv*dv;
      }
      vert_dist = sqrt(vert_dist);
      if (vert_dist>1.0) continue;

      std::cout << "matches neutrino vertex: vert_dist=" << vert_dist
      << " mctrack id=" << track.TrackID() << " pdg=" << track.PdgCode() << std::endl;

      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();
      std::vector<float> lepton_end(3);
      lepton_end[0] = last_step.X();
      lepton_end[1] = last_step.Y();
      lepton_end[2] = last_step.Z();
      std::cout << "lepton end=" << lepton_end[0] << "," << lepton_end[1] << "," << lepton_end[2] << std::endl;
      float norm = 0.;
      std::vector<float> lepton_dir(3);
      lepton_dir[0] = first_step.Momentum().Vect().X();
      lepton_dir[1] = first_step.Momentum().Vect().Y();
      lepton_dir[2] = first_step.Momentum().Vect().Z();
      for (int v=0; v<3; v++) norm += lepton_dir[v]*lepton_dir[v];
        norm = sqrt(norm);
      for (int v=0; v<3; v++) lepton_dir[v] /= norm;
        lepton_cosz = lepton_dir[2];
      lepton_phiz = atan2( lepton_dir[1], lepton_dir[0] );
      dwall_lepton = dwall( lepton_end, lepton_boundary );
    }//end of loop over mc tracks

    // Declare two vectors of 'larcv::Image2D'  that will be filled with empty images.
    std::vector < larcv::Image2D > init_markedimgs;
    std::vector < larcv::Image2D > final_markedimgs;

    // Test code....
    // Declare vectors for the cosmic pixels in the image
    // Taken from BMTrackCluster3D....
    for (size_t i = 0; i < 3; i++) {
      // I fill the image with the same dimensions as the image at this location in imgs_v.at(i) ('i' corresponds to the plane number).
      larcv::Image2D init_markedimg( imgs_v.at(i).meta() );
      init_markedimg.paint(0.0);
      init_markedimgs.emplace_back( std::move(init_markedimg) );
      final_markedimgs.emplace_back( std::move(init_markedimg) );

      // Cosmic Images 
    }
  }


    // count the pixels. determine if cosmic and neutrino are tagged. also if neutrino is in rois
    // we loop through the rows and cols
    int num_nupixels_inroi = 0;
    for (size_t p=0; p<3; p++) {
      for (size_t row=0; row<imgs_v.at(p).meta().rows(); row++) {
        for (size_t col=0; col<imgs_v.at(p).meta().cols(); col++) {
	  // check if this is a pixel of interest
          if ( imgs_v.at(p).pixel(row,col)<fthreshold ) continue;

          bool near_vertex = false;
	  // are we some radius from the vertex?
          if ( (int)row>=vertex_row-fvertex_radius && (int)row<=vertex_row+fvertex_radius
            && (int)col>=vertex_col[p]-fvertex_radius && (int)col<=vertex_col[p]+fvertex_radius ) {
            near_vertex = true;		
          }

          // above threshold. is it a neutrino pixel?
          const larcv::Image2D& segimg = segs_v.at(p);
          float x = imgs_v.at(p).meta().pos_x(col);
          float y = imgs_v.at(p).meta().pos_y(row);
          bool in_seg_image = false;
          int seg_row = -1;
          int seg_col = -1;
          if ( x>segs_v.at(p).meta().min_x() && x<segs_v.at(p).meta().max_x()
            && y>segs_v.at(p).meta().min_y() && y<segs_v.at(p).meta().max_y() ) {
            in_seg_image = true;
            seg_row = segs_v.at(p).meta().row(y);
            seg_col = segs_v.at(p).meta().col(x);
          }
          if ( in_seg_image && segs_v.at(p).pixel(seg_row,seg_col)>0 ) {

            nnu_pixels++;
            if (near_vertex)
              nvertex_pixels++;

	    // is it tagged?
            if ( stopmu_v.at(p).pixel(row,col)>0 || thrumu_v.at(p).pixel(row,col)>0 )  {
              nnu_tagged++;
              if ( near_vertex )
                nvertex_tagged++;
            }
            // is the neutrino pixel inside the ROI?
            for ( auto const& cand_roi : containedrois_v ) {
              float wired = imgs_v.at(p).meta().pos_x(col);
              float tick  = imgs_v.at(p).meta().pos_y(row);
              const larcv::ImageMeta& cand_roi_bb = cand_roi.BB().at(p);
              if ( cand_roi_bb.min_x()<wired && wired<cand_roi_bb.max_x() 
                && cand_roi_bb.min_y()<tick && tick<cand_roi_bb.max_y() )
                num_nupixels_inroi++;
            }
          }
          else {
	    // not a neutrino, so cosmic
            ncosmic_pixels++;
	    // is it tagged?
            if ( stopmu_v.at(p).pixel(row,col)>0 || thrumu_v.at(p).pixel(row,col)>0 ) 
              ncosmic_tagged++;
          
          }
        }//end of col loop
      }//end of row loop
    }//end of loop over planes for counting neutrino/cosmic pixels

    if ( nnu_pixels>0 )
      frac_inroi = float(num_nupixels_inroi)/float(nnu_pixels);
    else
      frac_inroi = 0.;
    std::cout << "fraction of neutrino pixels inside one of the rois: " << frac_inroi << std::endl;
    
    tree->Fill();

    if ( ientry>=100 )
      break;
  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
