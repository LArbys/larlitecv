#include "Base/DataCoordinator.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

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

  std::cout << "data[source] entries=" << dataco_source.get_nentries("larcv") << std::endl;
  std::cout << "data[thrumu] entries=" << dataco_thrumu.get_nentries("larcv") << std::endl;
  std::cout << "data[stopmu] entries=" << dataco_stopmu.get_nentries("larcv") << std::endl;

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
  float EnuGeV;        // neutrino energy in GeV
  float fdwall;        // dwall
  float dwall_lepton;  // dwll for end of lepton
  float frac_cosmic;   // fraction of cosmic tagged
  float frac_nu;       // fraction of neutrino pixels tagged
  float frac_vertex;   // fraction of vertex pixels tagged
  float primary_proton_ke; // ke of the proton driving from the hit nucleon
  float lepton_cosz;
  float lepton_phiz;
  float fpos[3];       // vertex
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
  tree->Branch("lepton_boundary",&lepton_boundary,"lepton_boundary/I");
  tree->Branch("vertex_boundary",&vertex_boundary,"vertex_boundary/I");
  tree->Branch("num_protons_over60mev", &num_protons_over60mev, "num_protons_over60mev/I");
  tree->Branch("EnuGeV",&EnuGeV,"EnuGeV/F");
  tree->Branch("dwall",&fdwall,"dwall/F");
  tree->Branch("dwall_lepton",&dwall_lepton,"dwall_lepton/F");
  tree->Branch("frac_cosmic",&frac_cosmic,"frac_cosmic/F");
  tree->Branch("frac_nu",&frac_nu,"frac_nu/F");
  tree->Branch("frac_vertex",&frac_vertex,"frac_vertex/F");
  tree->Branch("primary_proton_ke", &primary_proton_ke, "primary_proton_ke/F");
  tree->Branch("pos",fpos,"pos[3]/F");
  tree->Branch("lepton_cosz", &lepton_cosz, "lepton_cosz/F");
  tree->Branch("lepton_phiz", &lepton_phiz, "lepton_phiz/F");

  int nentries = dataco_stopmu.get_nentries("larcv");

  dataco_source.goto_entry(0,"larcv");
  dataco_thrumu.goto_entry(0,"larcv");
  dataco_stopmu.goto_entry(0,"larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco_stopmu.goto_entry(ientry,"larcv");

    dataco_stopmu.get_id(run,subrun,event);

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << std::endl;
      std::cout << " (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

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
    fpos[0] = fpos[1] = fpos[2] = 0.;
    num_protons_over60mev = 0;
    primary_proton_ke = 0.;

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
    }

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
