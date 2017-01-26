#ifndef __ENDPOINT_FILTER__
#define __ENDPOINT_FILTER__

#include <vector>
#include "BoundarySpacePoint.h"
#include "BoundaryEndPt.h"

// LArCV
#include "DataFormat/Image2D.h"
#include "dbscan/DBSCANAlgo.h"



namespace larlitecv {

  class EndPointFilter {
  public:
    EndPointFilter() {};
    virtual ~EndPointFilter() {};

		void filterEndPts( const std::vector< const BoundarySpacePoint* >& source, 
			const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes );

  	void removeBoundaryAndFlashDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
  		const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes );

  	void removeSameBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
	  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes );


  	void removeDiffBoundaryDuplicates( const std::vector< const BoundarySpacePoint* >& source, 
	  	const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<int>& endpoint_passes );

	  bool areEndPointsNearbyAndOnSameCluster( const larlitecv::BoundaryEndPt& pta, const larlitecv::BoundaryEndPt& ptb, 
	  	const larcv::Image2D& img, const larcv::Image2D& badch, const float radius_cm, dbscan::dbPoints& opt_cluster );



  };
  
}

#endif
