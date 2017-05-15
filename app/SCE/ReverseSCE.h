#ifndef __REVERSE_SCE_H__
#define __REVERSE_SCE_H__

#include <vector>
#include "TFile.h"
#include "TH3D.h"

namespace larlitecv {

  class ReverseSCE {
  public:
    ReverseSCE();
    virtual ~ReverseSCE();

    std::vector<double> getOriginalPos( const std::vector<double>& shiftedpos ) const;
    std::vector<float>  getOriginalPos( const std::vector<float>&  shiftedpos ) const;
    
  protected:

    TFile* rfile;
    TH3D* shift[3];

  };

}

#endif
