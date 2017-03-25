#include "extractTruthMethods.h"

#include <vector>
#include "TTree.h"

#include "DataFormat/mcpart.h"
#include "DataFormat/mctrajectory.h"

#include "dwall.h"

namespace larlitecv {

  void TruthData_t::clear() {
    // vertex info
    mode = 0;
    current = 0;
    vertex_boundary = 0;
    num_protons_over60mev = 0;
    fdwall = 0;
    EnuGeV = 0;
    for ( int i=0; i<3; i++ )
      pos[i] = 0;

    // lepton info
    num_protons_over60mev = 0;
    primary_proton_ke = 0.;
    lepton_boundary = 0;
    dwall_lepton = 0;
    lepton_cosz = 0;
    lepton_phiz = 0;

  }

  void TruthData_t::bindToTree( TTree* tree ) {

    // vertex truth
    tree->Branch("mode",&mode,"mode/I");
    tree->Branch("current",&current,"current/I");
    tree->Branch("vertex_boundary",&vertex_boundary,"vertex_boundary/I");
    tree->Branch("EnuGeV",&EnuGeV,"EnuGeV/F");
    tree->Branch("dwall",&fdwall,"dwall/F");
    tree->Branch("pos",pos,"pos[3]/F");

    // lepton truth
    tree->Branch("lepton_boundary",&lepton_boundary,"lepton_boundary/I");
    tree->Branch("num_protons_over60mev", &num_protons_over60mev, "num_protons_over60mev/I");
    tree->Branch("dwall_lepton",&dwall_lepton,"dwall_lepton/F");
    tree->Branch("lepton_cosz", &lepton_cosz, "lepton_cosz/F");
    tree->Branch("lepton_phiz", &lepton_phiz, "lepton_phiz/F");
    tree->Branch("primary_proton_ke", &primary_proton_ke, "primary_proton_ke/F");
  }

  void extractTruth( const larlite::event_mctruth& mctruth, const larlite::event_mctrack& mctrack, TruthData_t& data ) {
    // extract the truth quantities of interest
    extractVertexTruth( mctruth, data );
    extractLeptonTruth( mctruth, mctrack, data );
  }

  void extractVertexTruth( const larlite::event_mctruth& mctruth, TruthData_t& data ) {
    // extract the truth quantities of interest
    const larlite::mcnu& neutrino = mctruth.at(0).GetNeutrino();
    data.mode = neutrino.InteractionType();
    data.current = neutrino.CCNC();
    data.EnuGeV = neutrino.Nu().Momentum(0).E();
    const TLorentzVector& nu_pos = neutrino.Nu().Position();
    std::vector<float> fpos_v(3);
    std::vector<double> dpos(3);
    fpos_v[0] = nu_pos.X();
    fpos_v[1] = nu_pos.Y();
    fpos_v[2] = nu_pos.Z();
    data.fdwall = larlitecv::dwall(fpos_v, data.vertex_boundary);
    for (int i=0; i<3; i++)
      data.pos[i] = fpos_v[i];
  }

  void extractLeptonTruth( const larlite::event_mctruth& mctruth_v, const larlite::event_mctrack& mctrack_v, TruthData_t& data ) {

    // get the initial direction of the lepton
    const std::vector<larlite::mcpart>& particles = mctruth_v.at(0).GetParticles();
    bool found_lepton = false;
    int hit_nucleon_id = 2;
    int lepton_track_id = -1;
    int lepton_pid = -1;
    std::set<int> protonids;
    for ( auto  const& particle : particles ) {
      float KE = particle.Trajectory().front().E() - particle.Mass();
      if ( !found_lepton && (particle.PdgCode()==13 || particle.PdgCode()==-13 || particle.PdgCode()==11 || particle.PdgCode()==-11) ) {
        // found the lepton
        const larlite::mctrajectory& traj = particle.Trajectory();
        found_lepton = true;
        lepton_track_id = particle.TrackId();
        lepton_pid = particle.PdgCode();
        std::cout << "  lepton " << lepton_pid << ": E=" << traj.front().E() << " KE=" << KE << std::endl;
      }
      else if ( particle.PdgCode()==2212 ) {
        // look for the proton
        std::cout << "  a proton. p=" << particle.Momentum(0).Vect().Mag()
          << " E=" << particle.Trajectory().front().E() << " KE=" << KE
          << " status=" << particle.StatusCode()
          << " trackid=" << particle.TrackId() << " mother=" << particle.Mother()
          << std::endl;
        if ( particle.StatusCode()!=11 && KE>0.060 && protonids.find(particle.Mother())==protonids.end() )
          data.num_protons_over60mev++; // status 11 means from genie? threshold cut. check that it isn't from the original proton
        protonids.insert( particle.TrackId() );
      }
      else if ( particle.PdgCode()==14 || particle.PdgCode()==-14 || particle.PdgCode()==12 || particle.PdgCode()==-12 ) {
        std::cout << "  the neutrino (pdg=" << particle.PdgCode() << ") Enu=" << particle.Trajectory().front().E() << std::endl;
      }
      else {
        std::cout << "  pdg=" << particle.PdgCode()
          << " E=" << particle.Trajectory().front().E() << " KE=" << KE
          << " status=" << particle.StatusCode()
          << " end process=" << particle.EndProcess()
          << " trackid=" << particle.TrackId()
          << " mother=" << particle.Mother()
          << std::endl;
      }

      // stuff we are saving
      if ( (particle.PdgCode()==2212 || particle.PdgCode()==2112) && particle.StatusCode()==11 ) {
        hit_nucleon_id = particle.TrackId();
      }
      if ( particle.PdgCode()==2212 && particle.Mother()==hit_nucleon_id ) {
        data.primary_proton_ke = KE;
        protonids.insert(particle.TrackId());
      }

    }//end of particle track loop


    std::cout << "lepton track id = " << lepton_track_id << std::endl;
    std::cout << "num_protons_over60mev=" << data.num_protons_over60mev << std::endl;
    std::cout << "primary_proton_ke=" << data.primary_proton_ke << std::endl;

    // loop over MC tracks, find the neutrino lepton by matching vertex
    for ( auto const& track : mctrack_v ) {
      if ( std::abs(track.PdgCode())!=std::abs(lepton_pid)   ) continue;
      if ( track.size()==0 ) continue;
      const TLorentzVector& track_start = track.front().Position();
      std::vector<float> fstart(3);
      fstart[0] = track_start.X();
      fstart[1] = track_start.Y();
      fstart[2] = track_start.Z();

      float vert_dist = 0.;
      for (int v=0; v<3; v++) {
        float dv = data.pos[v]-fstart[v];
        vert_dist += dv*dv;
      }
      vert_dist = sqrt(vert_dist);
      std::cout << "  matches neutrino vertex: vert_dist=" << vert_dist << std::endl;

      const larlite::mcstep& first_step = track.front();
      const larlite::mcstep& last_step  = track.back();

      std::vector<float> lepton_end(3);
      lepton_end[0] = last_step.X();
      lepton_end[1] = last_step.Y();
      lepton_end[2] = last_step.Z();

      float norm = 0.;
      std::vector<float> lepton_dir(3);
      lepton_dir[0] = first_step.Momentum().Vect().X();
      lepton_dir[1] = first_step.Momentum().Vect().Y();
      lepton_dir[2] = first_step.Momentum().Vect().Z();
      for (int v=0; v<3; v++) norm += lepton_dir[v]*lepton_dir[v];
        norm = sqrt(norm);
      for (int v=0; v<3; v++) lepton_dir[v] /= norm;
        data.lepton_cosz = lepton_dir[2];
      data.lepton_phiz = atan2( lepton_dir[1], lepton_dir[0] );
      data.dwall_lepton = larlitecv::dwall( lepton_end, data.lepton_boundary );
    }//end of loop over mc tracks
  }


}