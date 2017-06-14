#include "ClusterGroupMatchingTypes.h"

namespace larlitecv {

  std::vector< larcv::Pixel2DCluster > ChargeVolume::GetPixelsInsideVolume( const std::vector<larcv::Image2D>& img_v,
    const std::vector<const ClusterGroup*>& clustergroups ) {
  	if ( img_v.size()!=m_clustergroups.size()){
  		throw std::runtime_error("number of image and plane clusters not the same");
  	}

  	std::vector< larcv::Pixel2DCluster > plane_pixels;
  	for ( size_t p=0; p<m_clustergroups.size(); p++) {

  		larcv::Pixel2DCluster pix_v;

      int num_pix_total = 0;
      int num_pix_invol = 0;

  		const ClusterGroup& group  = *(clustergroups.at(p));
  		const larcv::ImageMeta& meta = img_v.at(p).meta();

  		for ( auto const& slice : slices ) {
    		const WireInterval& winterval = slice.wire_intervals.at(p);
    		const std::vector<int>& rinterval = slice.row_interval;
    		for ( auto const& pixels : group.m_clusters_v ) {
    			for ( auto const& pix : pixels ) {
            num_pix_total++;
    				if ( (int)meta.pos_x(pix.X())>=winterval[0] && (int)meta.pos_x(pix.X())<=winterval[1]
    					   && (int)pix.Y()>=rinterval[0] && (int)pix.Y()<=rinterval[1] ) {
    					pix_v.push_back( pix );
              num_pix_invol++;
    				}
    			}
    		}
    	}
      std::cout << "fraction of pixels in volume: " << float(num_pix_invol)/float(num_pix_total) << " invol=" << num_pix_invol << " tot=" << num_pix_total << std::endl;
    	plane_pixels.emplace_back( std::move(pix_v) );
  	}

  	return plane_pixels;
  }

}

