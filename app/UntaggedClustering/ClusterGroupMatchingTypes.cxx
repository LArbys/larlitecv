#include "ClusterGroupMatchingTypes.h"

namespace larlitecv {

  std::vector< larcv::Pixel2DCluster > ChargeVolume::GetPixelsInsideVolume( const std::vector<larcv::Image2D>& img_v ) {
  	if ( img_v.size()!=m_clustergroups.size()){
  		throw std::runtime_error("number of image and plane clusters not the same");
  	}

  	std::vector< larcv::Pixel2DCluster > plane_pixels;
  	for ( size_t p=0; p<m_clustergroups.size(); p++) {

  		larcv::Pixel2DCluster pix_v;

  		const ClusterGroup& group  = *(m_clustergroups.at(p));
  		const larcv::ImageMeta& meta = img_v.at(p).meta();

  		for ( auto const& slice : slices ) {
    		const WireInterval& winterval = slice.wire_intervals.at(p);
    		const std::vector<int>& rinterval = slice.row_interval;
    		for ( auto const& pixels : group.m_clusters_v ) {
    			for ( auto const& pix : pixels ) {
    				if ( (int)meta.pos_x(pix.X())>=winterval[0] && (int)meta.pos_x(pix.X())<=winterval[1]
    					   && (int)pix.Y()>=rinterval[0] && (int)pix.Y()<=rinterval[1] ) {
    					pix_v.push_back( pix );
    				}
    			}
    		}
    	}
    	plane_pixels.emplace_back( std::move(pix_v) );
  	}

  	return plane_pixels;
  }

}

