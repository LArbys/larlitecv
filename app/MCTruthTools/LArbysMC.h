#ifndef __LARLITECV_LARBYSMC_H__
#define __LARLITECV_LARBYSMC_H__

#include <vector>
#include <string>

#include "TTree.h"

#include "DataFormat/storage_manager.h"

namespace larutil {
  class SpaceChargeMicroBooNE;
}

namespace larlitecv {

  /**
     \class ProcessBase
     User defined class LArbysImageMC ... these comments are used to generate
     doxygen documentation!
  */
  class LArbysMC  {

  public:

    /// Default constructor
    LArbysMC();

    /// Default destructor
    ~LArbysMC();

    //void configure(const PSet&);
    void initialize();
    bool process(larlite::storage_manager& mgr);
    void finalize();
    void Clear();

    void bindAnaVariables( TTree* );
    void printInteractionInfo();
    
  protected:

    std::string _name;
    
    std::string _producer_mctruth;
    std::string _producer_mctrack;
    std::string _producer_mcshower;    

    TTree* _mc_tree;

    larutil::SpaceChargeMicroBooNE* _psce;

    /// Event ID
    int _run;
    int _subrun;
    int _event;
    int _entry;

    std::string _rse_producer;

    // neutrino interaction info
    bool _neutrino_present;    
    int _current_type;
    int _interaction_type;
    int _genie_mode;

    float _Enu_true;
    float _vtx_t;
    float _vtx_x;
    float _vtx_y;
    float _vtx_z;
    float _vtx_sce_x;
    float _vtx_sce_y;
    float _vtx_sce_z;
    float _vtx_tick;
    float _vtx_wire[3];

    float _evis;
    float _evis_had;
    float _evis_vtx;
    float _evis_lep;

    int _hi_lep_pdg;
    float _hi_lep_e;
    
    /// Primary Particle Info
    int _parent_pdg;//primary particle pdg

    double _energy_deposit;
    double _energy_init;

    double _parent_x;
    double _parent_y;
    double _parent_z;
    double _parent_t;
    double _parent_px;
    double _parent_py;
    double _parent_pz;

    float _length_2d;

    int _nprimary;
    int _ntotal;

    int _nproton;
    int _nproton_60mev;    
    int _nlepton;
    int _nlepton_35mev;
    int _nmeson;
    int _nmeson_35mev;
    int _npi0;
    int _nshower;
    int _nneutron;
    int _1l1p0pi;

    double _dep_sum_lepton;
    double _ke_sum_lepton;

    double _dep_sum_proton;
    double _ke_sum_proton;

    double _dep_sum_meson;
    double _ke_sum_meson;

    double _dep_sum_shower;
    double _ke_sum_shower;

    double _dep_sum_neutron;
    double _ke_sum_neutron;

    /* geo2d::Vector<float> _start; //2d start point */
    /* geo2d::Vector<float> _dir;   //2d dir */

    std::vector<int>   _daughter_pdg_v;
    std::vector<int>   _daughter_trackid_v;
    std::vector<int>   _daughter_parenttrackid_v;

    std::vector<double> _daughter_energyinit_v;
    std::vector<double> _daughter_energydep_v;

    std::vector<std::vector<double> > _daughter_length_vv;
    std::vector<std::vector<double> > _daughter_2dstartx_vv;
    std::vector<std::vector<double> > _daughter_2dstarty_vv;
    std::vector<std::vector<double> > _daughter_2dendx_vv;
    std::vector<std::vector<double> > _daughter_2dendy_vv;
    std::vector<std::vector<double> > _daughter_2dcosangle_vv;

    std::vector<double> _daughterPx_v;
    std::vector<double> _daughterPy_v;
    std::vector<double> _daughterPz_v;
    std::vector<double> _daughterX_v;
    std::vector<double> _daughterY_v;
    std::vector<double> _daughterZ_v;
    std::vector<double> _daughter_length3d_v;

    double _true_proton_end_pt_x;
    double _true_proton_end_pt_y;
    double _true_proton_end_pt_z;

    double _proton_1st_pt_x;
    double _proton_1st_pt_y;
    double _proton_1st_pt_z;

    double _proton_last_pt_x;
    double _proton_last_pt_y;
    double _proton_last_pt_z;

    double _true_lepton_end_pt_x;
    double _true_lepton_end_pt_y;
    double _true_lepton_end_pt_z;

    double _lepton_1st_pt_x;
    double _lepton_1st_pt_y;
    double _lepton_1st_pt_z;

    double _lepton_last_pt_x;
    double _lepton_last_pt_y;
    double _lepton_last_pt_z;

  public:
    /// 2D Vertex Info
    std::vector<double> _vtx_2d_w_v;
    std::vector<double> _vtx_2d_t_v;


  protected:
    /// LARCV Image2D data
    /* std::vector<larcv::Image2D> _image_v; */
    /* ImageMeta _meta; */

    /// Filtering
    float _min_nu_dep_e;
    float _max_nu_dep_e;
    float _min_nu_init_e;
    float _max_nu_init_e;

    int _min_n_proton;
    int _min_n_neutron;
    int _min_n_lepton;
    int _min_n_meson;
    int _min_n_shower;

    float _min_proton_init_e;
    float _min_proton_dep;
    float _max_proton_dep;
    float _min_lepton_init_e;

    bool _check_vis;

    bool _do_not_reco;

    bool _selected;

    bool _select_signal;
    bool _select_background;

    struct aparticle{
      int trackid;
      int parenttrackid;
      float depeng;
      bool primary;

      bool daughterof (const aparticle& particle) const
      { return (this->parenttrackid == particle.trackid); }
    };


    struct entry_info{
      int run;
      int subrun;
      int event;
    };

    entry_info _entry_info;
    bool _is_signal;


  };

  /**
     \class larcv::LArbysMCFactory
     \brief A concrete factory class for larcv::LArbysMC
  */
  /* class LArbysMCProcessFactory : public ProcessFactoryBase { */
  /* public: */
  /*   /// ctor */
  /*   LArbysMCProcessFactory() { ProcessFactory::get().add_factory("LArbysMC",this); } */
  /*   /// dtor */
  /*   ~LArbysMCProcessFactory() {} */
  /*   /// creation method */
  /*   ProcessBase* create(const std::string instance_name) { return new LArbysMC(instance_name); } */

  /* }; */

}

#endif
/** @} */ // end of doxygen group
