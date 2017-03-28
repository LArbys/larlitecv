#include "FlashMatchMetricMethods.h"

#include <iostream>

namespace larlitecv {

  float CalculateFlashMatchChi2( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
    float& totpe_data, float& totpe_hypo, const float fudge_factor, const bool verbose ) {

    const int NPMTS = 32;
    float smallest_chi2 = -1.0;
    totpe_hypo = 0.;
    totpe_data = 0.;

    if ( verbose )
      std::cout << "CHI2 FLASH SCORE" << std::endl;
    for ( auto const& intime_flash : flash_data_v ) {
      float ll = 0.;
      float tot_pe = 0.;
      totpe_hypo = 0.;
      for (int i=0; i<NPMTS; i++) {
        float observed = intime_flash.PE(i);
        float expected = flash_hypothesis.PE(i)*fudge_factor;
        tot_pe += observed;
        totpe_hypo += expected;
        if ( observed>0 && expected==0 )
          expected = 1.0e-3;

        float dpe = (expected-observed);
        if ( observed>0 )
          dpe += observed*( log( observed ) - log( expected ) );
        ll += 2.0*dpe;
      }
      ll /= float(NPMTS);
      if ( smallest_chi2<0 || ll<smallest_chi2 ) {
        smallest_chi2 = ll;
        totpe_data = tot_pe;
      }

      if ( verbose ) {
        std::cout << "  [observed] ";
        for (int ich=0; ich<NPMTS; ich++ ) {
          std::cout << std::setw(5) << (int)intime_flash.PE(ich);
        }
        std::cout << " TOT=" << tot_pe << " LL=" << ll << std::endl;
      }
    }

    if ( verbose ) {
      std::cout << "  [expected] ";
      for ( int ich=0; ich<NPMTS; ich++ ) {
        std::cout << std::setw(5) << (int)(flash_hypothesis.PE(ich)*fudge_factor);
      }
      std::cout << " TOT[hypo]=" << totpe_hypo << " TOT[best data]=" << totpe_data << " BestLL=" << smallest_chi2 << std::endl;
    }

    return smallest_chi2;
  }

}
