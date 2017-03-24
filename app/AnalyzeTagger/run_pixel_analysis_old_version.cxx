#include "Base/DataCoordinator.h"

#include <iostream>
#include <string>
#include <sstream>

// ROOT
#include "TFile.h"
#include "TTree.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventROI.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larlitecv
#include "dwall.h"

// OpenCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"

int main( int nargs, char** argv ) {

  // PIXEL ANALYSIS

  if ( nargs!=2 ) {
    std::cout << "usage: ./run_pixel_analysis_old_version [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("AnalyzeTagger");

  enum { kSource=0, kCROIfile, kNumSourceTypes } SourceTypes_t;
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "AnalyzeTagger" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }

  // configuration parameters
  float fthreshold       = pset.get<float>("PixelThreshold");
  int tag_neighborhood   = pset.get<int>("TagNeighborhood",10);
  int fvertex_radius     = pset.get<int>("PixelRadius");
  int verbosity          = pset.get<int>("Verbosity",0);

  // setup output
  enum { kThruMu=0, kStopMu, kUntagged, kCROI, kNumStages } Stages_t;


  TFile* rfile = new TFile("output_pixel_analysis.root", "recreate");
  TTree* tree = new TTree("pixana", "Pixel-level analysis");

  // Event Index
  int run, subrun, event;

  // Truth Quantities about interaction
  int mode;                  // interaction mode
  int current;               // interaction cufrrent
  int lepton_boundary;       // index of boundary that lepton end point is nearest
  int vertex_boundary;       // index of boundary that vertex is nearest
  int num_protons_over60mev; // as named
  float EnuGeV;              // neutrino energy in GeV
  float fdwall;              // dwall
  float dwall_lepton;        // dwll for end of lepton
  float fpos[3];             // vertex
  float lepton_cosz;         // lepton dir
  float lepton_phiz;         // lepton dir
  float primary_proton_ke;   // ke of the proton driving from the hit nucleon

  // Truth Pixel Quantities
  int nthreshold_pixels[4]; // number of pixels above threshold
  int ncosmic_pixels[4];    // number of non-neutrino pixels
  int nnu_pixels[4];        // number of neutrino pixels
  int nvertex_pixels[4];    // number of neutrino pixels within some pixel radius of vertex

  // Tagged Pixel Quantities
  int ncosmic_tagged[kNumStages][4];       // number of non-neutrino pixels tagged
  int ncosmic_tagged_once[kNumStages][4];  // number of non-neutrino pixels tagged
  int ncosmic_tagged_many[kNumStages][4];  // number of non-neutrino pixels tagged
  int nnu_tagged[kNumStages][4];           // number of neutrino pixels tagged
  int nvertex_tagged[kNumStages][4];       // number of neutrino pixels within some pixel radius of vertex tagged
  std::stringstream s_arr;
  s_arr << "[" << (int)kNumStages << "][4]/I";

  // ROI quantities
  int num_rois;     // number of identified ROis
  int nnu_inroi[4]; // number of nu pixels contained in the CROI
  // float frac_cosmic;   // fraction of cosmic tagged
  // float frac_nu;       // fraction of neutrino pixels tagged
  // float frac_vertex;   // fraction of vertex pixels tagged
  // float frac_inroi;    // fraction of pixels contained in one of the rois

  // Event
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");

  // Truth
  tree->Branch("mode",&mode,"mode/I");
  tree->Branch("current",&current,"current/I");
  tree->Branch("lepton_boundary",&lepton_boundary,"lepton_boundary/I");
  tree->Branch("vertex_boundary",&vertex_boundary,"vertex_boundary/I");
  tree->Branch("num_protons_over60mev", &num_protons_over60mev, "num_protons_over60mev/I");
  tree->Branch("EnuGeV",&EnuGeV,"EnuGeV/F");
  tree->Branch("dwall",&fdwall,"dwall/F");
  tree->Branch("dwall_lepton",&dwall_lepton,"dwall_lepton/F");
  tree->Branch("pos",fpos,"pos[3]/F");
  tree->Branch("lepton_cosz", &lepton_cosz, "lepton_cosz/F");
  tree->Branch("lepton_phiz", &lepton_phiz, "lepton_phiz/F");
  tree->Branch("primary_proton_ke", &primary_proton_ke, "primary_proton_ke/F");

  tree->Branch("nthreshold_pixels", nthreshold_pixels, "nthreshold_pixels[4]/I");
  tree->Branch("ncosmic_pixels",    ncosmic_pixels,    "ncosmic_pixels[4]/I");
  tree->Branch("nnu_pixels",        nnu_pixels,        "nnu_pixels[4]/I");
  tree->Branch("nvertex_pixels",    nvertex_pixels,    "nvertex_pixels[4]/I");

  tree->Branch("ncosmic_tagged",      ncosmic_tagged,      std::string("ncosmic_tagged"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_once", ncosmic_tagged_once, std::string("ncosmic_tagged_once"+s_arr.str()).c_str() );
  tree->Branch("ncosmic_tagged_many", ncosmic_tagged_many, std::string("ncosmic_tagged_many"+s_arr.str()).c_str() );
  tree->Branch("nnu_tagged",          nnu_tagged,          std::string("nnu_tagged"+s_arr.str()).c_str() );
  tree->Branch("nvertex_tagged",      nvertex_tagged,      std::string("nvertex_tagged"+s_arr.str()).c_str() );

  tree->Branch("num_rois", &num_rois, "num_rois/I");
  tree->Branch("nnu_inroi", nnu_inroi, "nnu_inroi[4]" );
  //tree->Branch("frac_cosmic",&frac_cosmic,"frac_cosmic/F");
  //tree->Branch("frac_nu",&frac_nu,"frac_nu/F");
  //tree->Branch("frac_vertex",&frac_vertex,"frac_vertex/F");
  //tree->Branch("frac_inroi", &frac_inroi, "frac_inroi/F");


  int nentries = dataco[kCROIfile].get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(run,subrun,event);

    for ( int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
      if ( isrc!=kCROIfile ) {
        dataco[isrc].goto_event(run,subrun,event, "larcv");
      }
    }

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << ": (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    // initialize the output variables
    for (int p=0; p<4; p++) {
      ncosmic_pixels[p] = 0;
      nnu_pixels[0] = 0;
      nvertex_pixels[p] = 0;
      nthreshold_pixels[p] = 0;
      nnu_inroi[p] = 0;
      for (int istage=0; istage<kNumStages; istage++) {
	      ncosmic_tagged[istage][p] = 0;
	      ncosmic_tagged_once[istage][p] = 0;
	      ncosmic_tagged_many[istage][p] = 0;
	      nnu_tagged[istage][p] = 0;
	      nvertex_tagged[istage][p] = 0;
      }
    }
    mode = 0;
    current = 0;
    EnuGeV = 0.;
    fdwall = 0.;
    //frac_cosmic = 0.;
    //frac_vertex = 0.;
    //frac_nu = 0.;
    fpos[0] = fpos[1] = fpos[2] = 0.;
    num_protons_over60mev = 0;
    primary_proton_ke = 0.;

    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_segs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"modimg");
    larcv::EventImage2D* ev_badch  = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"gapchs");

    larcv::EventPixel2D* ev_pix[kNumStages] = {0};
    ev_pix[kThruMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"thrumupixels");
    ev_pix[kStopMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"stopmupixels");
    ev_pix[kUntagged]  = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"untaggedpixels");
    ev_pix[kCROI]      = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,"croipixels");

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>& segs_v   = ev_segs->Image2DArray();

    // get the result of the contained ROI analysis
    larcv::EventROI* ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    const std::vector<larcv::ROI>& containedrois_v = ev_contained_roi->ROIArray();
    num_rois = (int)containedrois_v.size();
    std::cout << "==ROIs==" << std::endl;
    for ( auto& roi : containedrois_v ) {
      std::cout << " roi: " << roi.dump();
    }

    // get other information, e.g. truth
    larlite::event_mctruth* ev_mctruth   = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
    larlite::event_mctrack* ev_mctrack   = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");

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
    fdwall = larlitecv::dwall(fpos_v, vertex_boundary);
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

    std::cout << "VERTEX (" << vertex_row << ", " << vertex_col[0] << "," << vertex_col[1] << "," << vertex_col[2] << ")" << std::endl;

    // get the initial direction of the lepton
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
	      std::cout << "  lepton E=" << traj.front().E() << " KE=" << KE << std::endl;
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
      dwall_lepton = larlitecv::dwall( lepton_end, lepton_boundary );
    }//end of loop over mc tracks

    // make thruth pixel counts
    // count the pixels. determine if cosmic and neutrino are tagged. also if neutrino is in rois
    // we loop through the rows and cols
    std::vector<larcv::Image2D> nupix_imgs_v;
    for (size_t p=0; p<3; p++) {
      // we create a neutrino pixel image, to make things easier downstream
      larcv::Image2D nupix_img( imgs_v.at(p).meta() );
      nupix_img.paint(0);
      for (size_t row=0; row<imgs_v.at(p).meta().rows(); row++) {
        for (size_t col=0; col<imgs_v.at(p).meta().cols(); col++) {
	        // check if this is a pixel of interest
          if ( imgs_v.at(p).pixel(row,col)<fthreshold ) continue;
	        nthreshold_pixels[p]++;
	        nthreshold_pixels[3]++;

          bool near_vertex = false;
          bool is_nu_pix = false;

	        // are we some radius from the vertex?
          if ( (int)row>=vertex_row-fvertex_radius && (int)row<=vertex_row+fvertex_radius
            && (int)col>=vertex_col[p]-fvertex_radius && (int)col<=vertex_col[p]+fvertex_radius ) {
            near_vertex = true;
            nvertex_pixels[p]++;
            nvertex_pixels[3]++;
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

            is_nu_pix = true;
            nnu_pixels[p]++;
            nnu_pixels[3]++;
            //if (near_vertex) {
            //  nvertex_pixels[p]++;
            //  nvertex_pixels[3]++;
	          //}

            // is the neutrino pixel inside the ROI?
            for ( auto const& cand_roi : containedrois_v ) {
              float wired = imgs_v.at(p).meta().pos_x(col);
              float tick  = imgs_v.at(p).meta().pos_y(row);
              const larcv::ImageMeta& cand_roi_bb = cand_roi.BB().at(p);
              if ( cand_roi_bb.min_x()<wired && wired<cand_roi_bb.max_x()
		            && cand_roi_bb.min_y()<tick && tick<cand_roi_bb.max_y() )
              {
                nnu_inroi[p]++;
                nnu_inroi[3]++;
	            }
            }
          }//if in semgment image
          else {
	          // not a neutrino, so cosmic
            ncosmic_pixels[p]++;
            ncosmic_pixels[3]++;
          }//end if cosmic

          if ( is_nu_pix )
            nupix_img.set_pixel( row, col, 1.0 );
          if ( near_vertex )
            nupix_img.set_pixel( row, col, 10.0 );

        }//end of col loop
      }//end of row loop
      nupix_imgs_v.emplace_back( std::move(nupix_img) );
    }//end of loop over planes for counting neutrino/cosmic pixels


    // Loop over track pixels
    // make left over image
    std::vector<larcv::Image2D> leftover_v;
    for ( auto const& img : imgs_v ) {
      larcv::Image2D leftover(img);
      leftover_v.emplace_back(leftover);
    }

    for ( int istage=0; istage<kNumStages; istage++ ) {
      for ( size_t p=0; p<3; p++ ) {
        // we need images to track the number of times pixels are Tagged
        larcv::Image2D img_counter( imgs_v.at(p).meta() );
        img_counter.paint(0);

        for ( auto const& pixcluster : ev_pix[istage]->Pixel2DClusterArray( p ) ) {

          std::set< std::vector<int> > tagged_set;
          for ( auto const& pix : pixcluster ) {

            for ( int dr=-tag_neighborhood; dr<=tag_neighborhood; dr++ ) {
              int r = pix.Y()+dr;
              if ( r<0 || r>=imgs_v.at(p).meta().rows() ) continue;
              for ( int dc=-tag_neighborhood; dc<=tag_neighborhood; dc++ ) {
                int c = pix.X()+dc;
                if ( c<0 || c>=imgs_v.at(p).meta().cols() ) continue;
                if ( imgs_v.at(p).pixel(r,c)<fthreshold ) continue;

                std::vector<int> pixel(2);
                pixel[0] = c;
                pixel[1] = r;
                tagged_set.insert( pixel );
              }
            }
          }//end of pixel loop
          //std::cout << "stage " << istage << " pix in cluster=" << pixcluster.size() << " tagged=" << tagged_set.size() << std::endl;
          for ( auto const& pix : tagged_set ) {
            float val = img_counter.pixel( pix[1], pix[0]) + 1.0;
            img_counter.set_pixel( pix[1], pix[0], val );
            leftover_v.at(p).set_pixel(pix[1],pix[0],0);
          }
        }//end of pix cluster loop

        // total the pixels
        for ( size_t r=0; r<img_counter.meta().rows(); r++) {
          for ( size_t c=0; c<img_counter.meta().cols(); c++ ) {

	          if ( imgs_v.at(p).pixel(r,c)<fthreshold ) continue;

            // is pixel cosmic or neutrino
            bool isnupix  = ( nupix_imgs_v.at(p).pixel(r,c)>=1.0  ) ? true : false;
            bool isvertex = ( nupix_imgs_v.at(p).pixel(r,c)>=10.0 ) ? true : false;

            int count = img_counter.pixel( r, c );

            if ( count>0 ) {
              if ( isnupix ) {
                nnu_tagged[istage][p]++;
                nnu_tagged[istage][3]++;
              }
              if ( isvertex ) {
                nvertex_tagged[istage][p]++;
                nvertex_tagged[istage][3]++;
              }
              if ( !isnupix && !isvertex ) {
                // is cosmic
                ncosmic_tagged[istage][p]++;
                ncosmic_tagged[istage][3]++;
                if (count==1)  {
                  ncosmic_tagged_once[istage][p]++;
                  ncosmic_tagged_once[istage][3]++;
                }
                else if ( count>1 ) {
                  ncosmic_tagged_many[istage][p]++;
                  ncosmic_tagged_many[istage][3]++;
                }
              }
            }
          }
        }
      }
    }

    /*
    int p = 0;
    for ( auto const& leftover : leftover_v ) {
      cv::Mat cvimg = larcv::as_mat_greyscale2bgr( leftover, fthreshold, 100.0 );
      std::stringstream ss;
      ss << "leftover_i" << ientry << "_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".jpg";
      cv::imwrite( ss.str(), cvimg );
      p++;
    }
    */

    tree->Fill();

  }//end of entry loop

  rfile->Write();
  return 0;

}//end of main
