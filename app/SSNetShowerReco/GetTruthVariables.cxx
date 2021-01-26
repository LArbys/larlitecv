#include "GetTruthVariables.h"

namespace larlitecv {
namespace ssnetshowerreco {

  /**
   * constructor.
   *
   * sets some defaults for parameters
   */
  GetTruthVariables::GetTruthVariables() {
    _mctruth_name = "generator";
    _mcshower_tree_name= "mcreco";
    _ssnet_shower_image_stem = "ubspurn";
    _thrumu_tree_name = "wire";
    clear();
  }

  void GetTruthVariables::initialize() {
    clear();
  }

  //write to outputfile
  void GetTruthVariables::finalize(){
    clear();
  }//end of finalize function

  void GetTruthVariables::clear() {
    _particle_PDG_v.clear();
    _particle_status_v.clear();
    _particle_momx_v.clear();
    _particle_momy_v.clear();
    _particle_momz_v.clear();
    _particle_E_v.clear();
    _gamma_energy_v.clear();

  }
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

  bool GetTruthVariables::process( larcv::IOManager& iolcv, larlite::storage_manager& ioll, int entry ) {

    // get truth information
    // store in root ana tree (if requested), store in json file
    // saving output particle information

    // clear result vectors first
    clear();

    //print out r,s,e (saved later)
    int _run    = ioll.run_id();
    int _subrun = ioll.subrun_id();
    int _event  = ioll.event_id();
    std::cout<<"(r,s,e)="<<_run<<","<<_subrun<<","<<_event<<std::endl;

    //load in inputs
    larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,  _mctruth_name );
    larlite::event_mcshower* ev_mcshower = ((larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,  _mcshower_tree_name));

    larcv::EventImage2D* ev_shower_score[3] = { nullptr };
    for ( size_t p=0; p<3; p++ ) {
      char treename[50];
      sprintf( treename, "%s_plane%d", _ssnet_shower_image_stem.c_str(), (int)p );
      ev_shower_score[p] =
        (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, treename );
    }
    larcv::EventImage2D* ev_thrumu
      = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, _thrumu_tree_name );
    const std::vector<larcv::Image2D>& thrumu_v = ev_thrumu->Image2DArray();

    int num = 0;
    //get particle info to save
    // first loop through showers for gamma info - used for shower reco test plot
    // for  (int shower =0;shower<ev_mcshower->size();shower++){
    //   double energy = ev_mcshower->at(shower).Start().E();
    //
    //   if (ev_mcshower->at(shower).PdgCode()==22){
    //     _gamma_energy_v.push_back(energy);
    //   }
    // }
    // std::cout<<_gamma_energy_v[0]<<" "<<_gamma_energy_v[0]<<std::endl;

    // loop through all final state paritcles
    for(int part =0;part<(int)ev_mctruth->at(0).GetParticles().size();part++){


      if (ev_mctruth->at(0).GetParticles().at(part).StatusCode() == 1){
        int pdg = ev_mctruth->at(0).GetParticles().at(part).PdgCode();
        _particle_PDG_v.push_back(pdg);
        num++;

        TLorentzVector momvec = ev_mctruth->at(0).GetParticles().at(part).Momentum();
        _particle_E_v.push_back(momvec.E());
        _particle_momx_v.push_back(momvec.Px());
        _particle_momy_v.push_back(momvec.Py());
        _particle_momz_v.push_back(momvec.Pz());

        // std::cout<<"["<<pdg<<","<<momvec.E()<<","<<momvec.Px()<<","<<momvec.Py()<<","<<momvec.Pz()<<"]\n";
      }//end of if status = 1 statement
    }//end of loop over particles
    _numparts = num;
    std::cout<<"total particles: "<<_numparts<<std::endl;

    // some code to look at ssnet average score to test

    // larcv::Image2D thru_y = thrumu_v[2];
    // float numpix = 0;
    // float totalval =0;
    // // now go through thru images, get average SSNet Score - just save yplane score
    // for (int r = 0; r<(int)thru_y.meta().rows();r++){
    //   for (int c= 0; c<(int)thru_y.meta().cols();c++){
    //   // check if pass cosmic tagger
    //     if (true){
    //       float ssnetshowerscore = ev_shower_score[2]->Image2DArray()[0].pixel(r,c);
    //       numpix+=1.0;
    //       totalval+=ssnetshowerscore;
    //     }
    //   }
    // }
    // _averageSSNetShowerScore = totalval/numpix;
    // std::cout<<"Average SSNet Score: "<<_averageSSNetShowerScore<<std::endl;

    return true;
  }//end of process function

}
}
