#include "Base/DataCoordinator.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/mctruth.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "TFile.h"
#include "TTree.h"

float dwall( const std::vector<float>& pos ) {

  float dx1 = fabs(pos[0]);
  float dx2 = fabs(260.0-pos[0]);
  float dy1 = fabs(118.0-pos[1]);
  float dy2 = fabs(-118.0-pos[1]);
  float dz1 = fabs(pos[2]);
  float dz2 = fabs(1036.0-pos[2]);

  float dwall = 1.0e9;
  if ( dx1<dwall ) dwall = dx1;
  if ( dx2<dwall ) dwall = dx2;
  if ( dy1<dwall ) dwall = dy1;
  if ( dy2<dwall ) dwall = dy2;
  if ( dz1<dwall ) dwall = dz1;
  if ( dz2<dwall ) dwall = dz2;

  return dwall;

}

int main( int nargs, char** argv ) {

  // run pixel analysis. use 

	// we need to have a data coordinator for each stage because the number of entries could be different.
	// we'll coordinate by using event,subrun,run information
	std::string data_folder = "~/data/larbys/cosmic_tagger/mcc7_bnbcosmic/";

  larlitecv::DataCoordinator dataco_source;
  dataco_source.add_inputfile( data_folder+"/output_larcv.root", "larcv" ); // segment image/original image
  dataco_source.add_inputfile( data_folder+"/output_larlite.root", "larlite"); //source larlite file

  larlitecv::DataCoordinator dataco_thrumu;
	dataco_thrumu.add_inputfile( data_folder+"/output_larcv_testmcbnbcosmic_signalnumu.root", "larcv" ); // thrumu-tagger info, larcv
  dataco_thrumu.add_inputfile( data_folder+"/output_larlite_testmcbnbcosmic_signalnumu.root", "larlite"); //thrumu-tagged info, larlite

  larlitecv::DataCoordinator dataco_stopmu;
  dataco_stopmu.add_inputfile( data_folder+"/output_stopmu_larcv_p1.root", "larcv" ); //stopmu-tagger output


  dataco_source.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_thrumu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );
  dataco_stopmu.configure( "config.cfg", "StorageManager", "IOManager", "PixelAnalysis" );

  dataco_source.initialize();
  dataco_thrumu.initialize();
  dataco_stopmu.initialize();

  // configuration parameters
  larcv::PSet cfg = larcv::CreatePSetFromFile( "config.cfg" );
  larcv::PSet pixana_cfg = cfg.get<larcv::PSet>("PixelAnalysis");
  float fthreshold   = pixana_cfg.get<float>("PixelThreshold");
  int fvertex_radius =  pixana_cfg.get<int>("PixelRadius");

  // setup output

  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
  TTree* tree = new TTree("pixana", "Pixel-level analysis");
  int ncosmic_pixels; // number of non-neutrino pixels
  int nnu_pixels;     // number of neutrino pixels
  int nvertex_pixels; // number of neutrino pixels within some pixel radius of vertex
  int ncosmic_tagged; // number of non-neutrino pixels tagged
  int nnu_tagged;     // number of neutrino pixels tagged
  int nvertex_tagged; // number of neutrino pixels within some pixel radius of vertex tagged
  int mode;           // interaction mode
  int current;        // interaction cufrrent
  float EnuGeV;       // neutrino energy in GeV
  float fdwall;        // dwall
  float frac_cosmic;  // fraction of cosmic tagged
  float frac_nu;      // fraction of neutrino pixels tagged
  float frac_vertex;  // fraction of vertex pixels tagged
  int run, subrun, event;
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");
	tree->Branch("ncosmic_pixels",&ncosmic_pixels,"ncosmic_pixels/I");
	tree->Branch("nnu_pixels",&nnu_pixels,"nnu_pixels/I");
	tree->Branch("nvertex_pixels",&nvertex_pixels,"nvertex_pixels/I");
	tree->Branch("ncosmic_tagged",&ncosmic_tagged,"ncosmic_tagged/I");
	tree->Branch("nnu_tagged",&nnu_tagged,"nnu_tagged/I");
	tree->Branch("nvertex_tagged",&nvertex_tagged,"nvertex_tagged/I");
	tree->Branch("mode",&mode,"mode/I");
	tree->Branch("current",&current,"current/I");
	tree->Branch("EnuGeV",&EnuGeV,"EnuGeV/F");
	tree->Branch("dwall",&fdwall,"dwall/F");
	tree->Branch("frac_cosmic",&frac_cosmic,"frac_cosmic/F");
  tree->Branch("frac_nu",&frac_nu,"frac_nu/F");
  tree->Branch("frac_vertex",&frac_vertex,"frac_vertex/F");

  int nentries = dataco_stopmu.get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

  	dataco_stopmu.goto_entry(ientry,"larcv");

  	dataco_stopmu.get_id(run,subrun,event);

  	std::cout << "entry " << ientry << std::endl;
  	std::cout << " (r,s,e)=(" << run << ", " << " , " << subrun << ", " << event << ")" << std::endl;

  	dataco_thrumu.goto_event(run,subrun,event,"larcv");
  	dataco_source.goto_event(run,subrun,event,"larcv");

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

  	// ok now to do damage

  	// get the original, segmentation, and tagged images
  	larcv::EventImage2D* ev_segs   = (larcv::EventImage2D*)dataco_source.get_larcv_data(larcv::kProductImage2D,"segment");
  	larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"modimgs");
  	//larcv::EventImage2D* ev_badch  = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"badch");	
  	larcv::EventImage2D* ev_thrumu = (larcv::EventImage2D*)dataco_thrumu.get_larcv_data(larcv::kProductImage2D,"marked3d");
  	larcv::EventImage2D* ev_stopmu = (larcv::EventImage2D*)dataco_stopmu.get_larcv_data(larcv::kProductImage2D,"stopmu");  	

		const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
		//const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
		const std::vector<larcv::Image2D>& segs_v   = ev_segs->Image2DArray();
		const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();
		const std::vector<larcv::Image2D>& stopmu_v = ev_stopmu->Image2DArray();

		// get other information, e.g. truth
		larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)dataco_source.get_larlite_data(larlite::data::kMCTruth,"generator");

		// extract the truth quantities of interest
    const larlite::mcnu& neutrino = ev_mctruth->at(0).GetNeutrino();

	  mode = neutrino.InteractionType();
	  current = neutrino.CCNC();
	  EnuGeV = neutrino.Nu().Momentum(0).E();
    const TLorentzVector& nu_pos = neutrino.Nu().Position();
    std::vector<float> fpos(3);
    std::vector<double> dpos(3);
    fpos[0] = nu_pos.X();
    fpos[1] = nu_pos.Y();
    fpos[2] = nu_pos.Z();
    dpos[0] = nu_pos.X();
    dpos[1] = nu_pos.Y();
    dpos[2] = nu_pos.Z();
    fdwall = dwall(fpos);

    // get the vertex in the pixel coordinates
    std::vector<int> wid(3,-1);
    std::vector<int> vertex_col(3,-1);
		for (size_t p=0; p<3; p++) {
			wid[p] = ::larutil::Geometry::GetME()->WireCoordinate( dpos, p );
			if ( wid[p]>=0 && wid[p]<3456 )
				vertex_col[p] = imgs_v.at(p).meta().col(wid[p]);
		}
		float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
		float ftick = nu_pos[0]/cm_per_tick + 3200.0;
		int vertex_row = -1;
		if ( ftick >= imgs_v.at(0).meta().min_y() && ftick<=imgs_v.at(0).meta().max_y() )
			vertex_row = imgs_v.at(0).meta().row( ftick );

		// count the pixels. we loop through the rows and cols
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
						}
						else {
							// not a neutrino, so cosmic
							ncosmic_pixels++;
							// is it tagged?
							if ( stopmu_v.at(p).pixel(row,col)>0 || thrumu_v.at(p).pixel(row,col)>0 ) 
								ncosmic_tagged++;
						}
					}
				}
			}
			tree->Fill();

			if ( ientry>=100 )
				break;
		}

	rfile->Write();
  return 0;
};
