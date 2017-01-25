#include "DBClusterMergerAlgo.h"

#include "AStarDirAlgo.h"

namespace larlitecv {

	bool DBClusterMergerAlgo::willClustersMergeThroughBadchs( const int clusterid_a, const int clusterid_b, 
			const dbscan::dbPoints& pixels, const dbscan::dbscanOutput& cluster_info, 
      const larcv::Image2D& img, const larcv::Image2D& badch_img, const DBClusterMergerPars_t& pars ) {
		// we use this function when we want to see if two clusters (produced by dbscan::DBSCANAlgo) can be connect by badchs
		// we use the A* algorithm to try and draw a path from one end of the one cluster to the other end of the other cluster.
		// the A* algo has the ability to path-find across badch regions. It also marks the jump as such as well.
		// kind of expensive routine

		// first we scan each cluster and find it's extrema
		dbscan::ClusterExtrema extrema_a = dbscan::ClusterExtrema::FindClusterExtrema( clusterid_a, cluster_info, pixels );
		dbscan::ClusterExtrema extrema_b = dbscan::ClusterExtrema::FindClusterExtrema( clusterid_b, cluster_info, pixels );

		// becuz I am lazy
		typedef dbscan::ClusterExtrema::Extrema_t ex_t;

		// we identify the end points that are the furtherest away from each other
		int bestcombo[2] = {0,0};
		float max_dist = -1;
		for (int ia=0; ia<dbscan::ClusterExtrema::kNumExtrema; ia++ ) {
			for (int ib=0; ib<dbscan::ClusterExtrema::kNumExtrema; ib++ ) {			

				float dc = extrema_a.extrema((ex_t)ia)[0] = extrema_b.extrema((ex_t)ib)[0];
				float dr = extrema_a.extrema((ex_t)ia)[1] = extrema_b.extrema((ex_t)ib)[1];

				float dist = 0.;
				dist += sqrt( dc*dc + dr*dr );

				if ( max_dist<0 || max_dist < dist ) {
					max_dist = dist;
					bestcombo[0] = ia;
					bestcombo[1] = ib;
				}
			}
		}

		// could not find best combo?
		if ( max_dist==-1 )
			return false;

		AStarDirAlgoConfig cfg;
		cfg.astar_threshold.resize( 3, pars.astar_pixel_threshold );
		cfg.astar_neighborhood.resize( 3, pars.astar_neighborhood );
		cfg.astar_start_padding = pars.astar_start_pad;
		cfg.astar_end_padding   = pars.astar_end_pad;		
		cfg.image_padding       = pars.astar_image_padding;

		AStarDirAlgo algo( cfg );
		int start_col = extrema_a.extrema((ex_t)bestcombo[0])[0];
		int start_row = extrema_a.extrema((ex_t)bestcombo[0])[1];		
		int end_col   = extrema_b.extrema((ex_t)bestcombo[1])[0];
		int end_row   = extrema_b.extrema((ex_t)bestcombo[1])[1];		
		algo.setBadChImage( badch_img );
		std::vector<AStarDirNode> path = algo.findpath( img, start_row, start_col, end_row, end_col, pars.astar_pixel_threshold, false );

		if ( path.back().row==end_row && path.back().col==end_col ) {
			// success!
			return true;
		}

		return false;
	}

}// end of namespace