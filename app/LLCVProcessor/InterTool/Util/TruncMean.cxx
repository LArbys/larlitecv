#ifndef TRUNCMEAN_CXX
#define TRUNCMEAN_CXX

#include "TruncMean.h"
#include <cassert>
#include <stdexcept>
#include <numeric>
namespace llcv {

  float TruncMean::CalcIterativeTruncMean(std::vector<float> v, const size_t& nmin,
					  const size_t& nmax, const size_t& currentiteration,
					  const size_t& lmin,
					  const float& convergencelimit,
					  const float& nsigma, const float& oldmed)
  {

    auto const& mean = Mean(v);
    auto const& med  = Median(v);
    auto const& rms  = RMS(v);

    // if the vector length is below the lower limit -> return
    if (v.size() < lmin)
      return mean;

    // if we have passed the maximum number of iterations -> return
    if (currentiteration >= nmax)
      return mean;

    // if we passed the minimum number of iterations and the mean is close enough to the old value
    float fracdiff = fabs(med-oldmed) / oldmed;
    if ( (currentiteration >= nmin) && (fracdiff < convergencelimit) )
      return mean;
  
    // if reached here it means we have to go on for another iteration

    // cutoff tails of distribution surrounding the mean
    // use erase-remove : https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
    // https://stackoverflow.com/questions/17270837/stdvector-removing-elements-which-fulfill-some-conditions
    v.erase( std::remove_if( v.begin(), v.end(), 
			     [med,nsigma,rms](const float& x) { return ( (x < (med-nsigma*rms)) || (x > (med+nsigma*rms)) ); }), // lamdda condition for events to be removed
	     v.end());
  
    return CalcIterativeTruncMean(v, nmin, nmax, lmin, currentiteration+1, convergencelimit, nsigma, med);
  }

  void TruncMean::CalcTruncMeanProfile(const std::vector<float>& rr_v, const std::vector<float>& dq_v,
				       std::vector<float>& dq_trunc_v, const float& nsigma)
  {

    // how many points to sample 
    int Nneighbor = (int)(_rad * 3 * 2);

    dq_trunc_v.clear();
    dq_trunc_v.reserve( rr_v.size() );

    int Nmax = dq_v.size() - 1;

    for (size_t n=0; n < dq_v.size(); n++) {

      // current residual range
      float rr = rr_v.at(n);

      int nmin = n - Nneighbor;
      int nmax = n + Nneighbor;

      if (nmin < 0) nmin = 0;
      if (nmax > Nmax) nmax = Nmax;

      // vector for local dq values
      std::vector<float> dq_local_v;

      for (int i=nmin; i < nmax; i++) {
      
	float dr = rr - rr_v[i];
	if (dr < 0) dr *= -1;

	if (dr > _rad) continue;

	dq_local_v.push_back( dq_v[i] );
      
      }// for all ticks we want to scan

      if (dq_local_v.size() == 0) {
	dq_trunc_v.push_back( dq_v.at(n) );
	continue;
      }
    
      // calculate median and rms
      float median = Median(dq_local_v);
      float rms    = RMS(dq_local_v);

      float truncated_dq = 0.;
      int npts = 0;
      for (auto const& dq : dq_local_v) {
	if ( ( dq < (median+rms * nsigma) ) && ( dq > (median-rms * nsigma) ) ){
	  truncated_dq += dq;
	  npts += 1;
	}
      }

      float ratio = truncated_dq;
	
      if (npts != 0)
	ratio /= (float) npts;

      dq_trunc_v.push_back(ratio);
    }// for all values

    return;
  }

  float TruncMean::Mean(const std::vector<float>& v)
  {
    if (v.empty())
      throw std::runtime_error("Empty input to Median");

    if (v.size() == 1) 
      return v.front();

    float mean = 0.;
    for (auto const& n : v) mean += n;
    mean /= (float) v.size();
  
    return mean;
  }

  float TruncMean::Median(const std::vector<float>& v)
  {
    if (v.empty())
      throw std::runtime_error("Empty input to Median");

    if (v.size() == 1) 
      return v.front();
  
    std::vector<float> vcpy = v;

    std::sort(vcpy.begin(), vcpy.end());

    float index = (float)(vcpy.size()) / 2.0;
    int idx = (size_t) index;

    if ((idx < 0) or (idx >= (int)vcpy.size()))
      throw std::runtime_error("Median index is outside of v range");
    
    float median = vcpy[(size_t)index];

    return median;
  }

  float TruncMean::RMS(const std::vector<float>& v)
  {

    if (v.empty())
      throw std::runtime_error("Empty input to RMS");

    if (v.size()==1) 
      return 0.0;
    
    float avg = 0.;
    for (auto const& val : v) avg += val;

    avg /= (float) v.size();

    float rms = 0.;

    for (auto const& val : v) rms += (val-avg)*(val-avg);

    float denom = (float)v.size() - 1;

    rms = sqrt( rms / denom );

    return rms;
  }
  
  float TruncMean::Slope(const std::vector<float>& x,
			 const std::vector<float>& y)
  {
    assert (x.size() == y.size());

    float n     = x.size();
    float s_x   = std::accumulate(x.begin(), x.end(), 0.0);
    float s_y   = std::accumulate(y.begin(), y.end(), 0.0);
    float s_xx  = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    float s_xy  = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    float slope = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

    return slope;
  }

  void TruncMean::Linearize(const std::vector<float>& v1, std::vector<float>& v2) {
    assert (v1.size() == v2.size());

    size_t t  = 0;
    size_t ts = 0;

    while(t < v2.size()) {
      
      if (v1[t] != 0)  {
	v2[t] = v1[t];

	if (ts != 0) {
	  size_t step  = ts;
	  size_t start = ts-1;
	  while(step <= t) {
	    v2[step]  = (v2[t] - v2[start]);
	    v2[step] /= (float)(t - start);
	    v2[step] *= (float)(step - start);
	    v2[step] += v2[start];
	    ++step;
	  }
	}
	ts = 0;
      }
      else { 
	if (ts == 0) 
	  ts = t;
      }
      ++t;
    }

    return;
  }
}
#endif
