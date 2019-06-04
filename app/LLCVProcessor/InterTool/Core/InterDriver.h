#ifndef __INTERDRIVER_H__
#define __INTERDRIVER_H__

#include "InterAnaBase.h"
#include "InterSelBase.h"

#include "InterTTreeManager.h"
#include "InterDataManager.h"
#include "InterImageManager.h"

#include "InterToolTypes.h"

#include "TFile.h"

namespace llcv {

  class InterModule;
  class InterMichel;
  class InterPMT;
  
  class InterDriver : public llcv_base {

    friend class InterModule;
    friend class InterMichel;
    friend class InterPMT;

  public:
  
  InterDriver(std::string name="InterDriver")
    :   llcv_base(name)
      ,_fout(nullptr)
    {}

    ~InterDriver() {}

  private:

    TFile* _fout;

    std::string _fout_fname;

    int _run;
    int _subrun;
    int _event;
    int _vtxid;
    
    // Next() for next vertex
    InterTTreeManager _tree_mgr;

    // data manager per vertex
    std::vector<InterDataManager> _data_mgr_v;

    // image manager per event
    InterImageManager _img_mgr;

    // analysis to be run per vertex or per event
    std::vector<InterSelBase*> _sel_base_v;
    
  public:

    //
    // Public configuration
    //
    void AddSelection(InterSelBase* sbase);
    void SetOutputFilename(std::string fout_fname);
    void AttachInterFile(const std::string& fname,const std::string& tname);
    
  private:

    //
    // InterModule directives
    //
    void Configure(const larcv::PSet& cfg);
    void Initialize();
    void Process();
    void Finalize();

    //
    // InterModule attachments (can change)
    //
    void AttachImage(const std::vector<larcv::Image2D>& img_v, InterImageType itype);
    
    size_t AttachVertex   (const larlite::vertex* vertex);
    size_t AttachPGraph   (size_t vtxid, const larcv::PGraph* pgraph);
    size_t AttachParticles(size_t vtxid, const larcv::PGraph* pgraph, const larcv::EventPixel2D* ev_pix);
    size_t AttachOpFlash  (size_t vtxid, const larlite::opflash* opflash);
    size_t AttachTrack    (size_t vtxid, const larlite::track* track);
    size_t AttachShower   (size_t vtxid, const larlite::shower* shower);
    size_t AttachCluster  (size_t vtxid, size_t shrid, const larlite::cluster* cluster);
    size_t AttachHit      (size_t vtxid, size_t cluid, const larlite::hit* hit);
    size_t AttachHit      (size_t vtxid, const larlite::hit* hit);
    size_t AttachOpDigit  (size_t vtxid, const larlite::opdetwaveform* opdigit);

    void Reset();
    void Dump();

  };
  
  
}

#endif
