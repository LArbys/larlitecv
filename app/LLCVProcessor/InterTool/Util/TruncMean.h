/**
 * \file TruncMean.h
 *
 * \ingroup 3DMichel
 * 
 * \brief Class def header for a class TruncMean
 *
 * @author david caratelli [davidc@fnal.gov]
 * Written 08/02/2018.
 */

#ifndef TRUNCMEAN_H
#define TRUNCMEAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>

/**
   \class TruncMean
   The truncated mean class allows to compute the following quantities
   1) the truncated mean profile of an ordered vector of values, such as 
   the charge profile along a particle's track.
   To create such a profile use the function CalcTruncMeanProfile()
   2) Get the truncated mean value of a distribution. This function
   iteratively hones in on the truncated mean of a distribution by
   updating the mean and cutting the tails after each iteration.
   For this functionality use CalcIterativeTruncMean()
   doxygen documentation!
*/

static const float kINVALID_FLOAT = std::numeric_limits<float>::max();
namespace llcv {

  class TruncMean{

  public:

    /// Default constructor
    TruncMean(){}

    /// Default destructor
    ~TruncMean(){}

    /**
       @brief Given residual range and dq vectors return truncated local dq.
       Input vectors are assumed to be match pair-wise (nth entry in rr_v
       corresponds to nth entry in dq_v vector).
       Input rr_v values are also assumed to be ordered: monotonically increasing
       or decreasing.
       For every dq value a truncated linear dq value is calculated as follows:
       0) all dq values within a rr range set by the class variable _rad are selected.
       1) the median and rms of these values is calculated.
       2) the subset of local dq values within the range [median-rms, median+rms] is selected.
       3) the resulting local truncated dq is the average of this truncated subset.
       @input std::vector<float> rr_v -> vector of x-axis coordinates (i.e. position for track profile)
       @input std::vector<float> dq_v -> vector of measured values for which truncated profile is requested
       (i.e. charge profile of a track)
       @input std::vector<float> dq_trunc_v -> passed by reference -> output stored here
       @input float nsigma -> optional parameter, number of sigma to keep around RMS for TM calculation
    */
    void CalcTruncMeanProfile(const std::vector<float>& rr_v, const std::vector<float>& dq_v,
			      std::vector<float>& dq_trunc_v, const float& nsigma = 1);

    /**
       @brief Iteratively calculate the truncated mean of a distribution
       @input std::vector<float> v -> vector of values for which truncated mean is asked
       @input size_t nmin -> minimum number of iterations to converge on truncated mean
       @input size_t nmax -> maximum number of iterations to converge on truncated mean
       @input size_t lmin -> minimum number of entries in vector before exiting and returning current value
       @input size_t currentiteration -> current iteration
       @input float convergencelimit -> fractional difference between successive iterations
       under which the iteration is completed, provided nmin iterations have occurred.
       @input nsigma -> number of sigma around the median value to keep when the distribution is trimmed.
    */
    float CalcIterativeTruncMean(std::vector<float> v, const size_t& nmin,
				 const size_t& nmax, const size_t& currentiteration,
				 const size_t& lmin,
				 const float& convergencelimit,
				 const float& nsigma, const float& oldmed = kINVALID_FLOAT);

    /**
       @brief Set the smearing radius over which to take hits for truncated mean computaton.
    */
    void setRadius(const float& rad) { _rad = rad; }

    //
    // some shits from vic
    //
    void Linearize(const std::vector<float>& v1, std::vector<float>& v2);

    float Mean  (const std::vector<float>& v);
    float Median(const std::vector<float>& v);
    float RMS   (const std::vector<float>& v);
    float Slope (const std::vector<float>& x,
		 const std::vector<float>& y);

  private:
    /**
       Smearing radius over which charge from neighboring hits is scanned to calculate local
       truncated mean
    */
    double _rad;


  };
}
#endif
/** @} */ // end of doxygen group 
