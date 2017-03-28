#include "FlashMatchMetricMethods.h"

#include <iostream>
#include <cmath>

#include "LArUtil/Geometry.h"

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

  // ====================================================================================================
  // Shape only 2D Unbinned Log-likelihood

  float CalculateShapeOnlyUnbinnedLL( const std::vector<larlite::opflash>& flash_data_v, const larlite::opflash& flash_hypothesis,
    float& totpe_data, float& totpe_hypo, const bool verbose ) {
    // here we use the flash hypo to make a 2D gassian likelihood. Simply use mean and covariance to get shape in (y,z)
    // then we
    float mean[2] = {0.0};
    float tot_weight = 0;
    totpe_data = -1.0;
    totpe_hypo = -1.0;

    std::vector<double> xyz(3,0.0);
    for (int ich=0; ich<32; ich++ ) {
      larutil::Geometry::GetME()->GetOpChannelPosition( ich, xyz );
      mean[0] += xyz[1]*flash_hypothesis.PE(ich);
      mean[1] += xyz[2]*flash_hypothesis.PE(ich);
      tot_weight = flash_hypothesis.PE(ich);
    }
    for (int i=0; i<2; i++)
      mean[i] /= tot_weight;

    float cov[2][2] = {0};
    for (int ich=0; ich<32; ich++) {
      larutil::Geometry::GetME()->GetOpChannelPosition( ich, xyz );
      float dx[2];
      dx[0] = xyz[1] - mean[0];
      dx[1] = xyz[2] - mean[1];
      for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++ ) {
          cov[i][j] += dx[i]*dx[j]*flash_hypothesis.PE(ich)/tot_weight;
        }
      }
    }
    totpe_hypo = tot_weight;

    float det = cov[0][0]*cov[1][1] - cov[1][0]*cov[0][1];
    float invcov[2][2] = { { cov[1][1]/det, -cov[0][1]/det},
                           {-cov[1][0]/det,  cov[0][0]/det } };
    float normfactor = 2*3.141159*det;

    float min_neglnl = 0.0;
    for ( int iflash=0; iflash<(int)flash_data_v.size(); iflash++ ) {
      const larlite::opflash& dataflash = flash_data_v.at(iflash);
      float neglnl = 32*0.5*log(normfactor);
      float totweight = 0.;
      for ( int ich=0; ich<32; ich++) {
        larutil::Geometry::GetME()->GetOpChannelPosition( ich, xyz );
        float mahadist = 0.;
        float yz[2] = { xyz[1], xyz[2] };
        for (int i=0; i<2; i++ ) {
          for (int j=0; j<2; j++) {
            mahadist += (yz[i]-mean[i])*invcov[i][j]*(yz[j]-mean[j]);
          }
        }
        neglnl += 0.5*mahadist*dataflash.PE(ich);
        totweight += dataflash.PE(ich);
      }
      neglnl /= totweight;
      if ( iflash==0 || neglnl < min_neglnl ) {
        min_neglnl = neglnl;
        totpe_data = totweight;
      }
    }

    return min_neglnl;
  }



}
