#include "ContainedROI.h"

namespace larlitecv {

	ContainedROI::ContainedROI( const ContainedROIConfig& config ) {
		m_config = config;
	} 


 	std::vector<larcv::ROI> ContainedROI::SelectROIs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v,
    	const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v, 
    	const larlite::event_opflash& opflash_v ) {

 		// goal is to return ROIs around contained clusters that will serve as potential neutrino candidates

 		// path
 		// cluster untagged pixels
 		// on each cluster, determine extrema points
 		// use extrema points to try to match up such clusters across three planes: look for same time, consistent 3D positions, and proportionate charge
 		// next, use z-position of trigger flash to collect all consistent thrumu, stopmu, and untagged clusters
 		// for stopmu clusters with flash ends, we check if trigger flash more consistent than tagging flash. if so, we convert it into an untagged cluster
 		// for all untagged clusters that match the trigger flash by some threshold, this gets packaged as an ROI and returned

 		// collect untagged pixels into its own image
 		std::vector<larcv::Image2D> untagged_v;
 		for ( size_t p=0; p<img_v.size(); p++ ) {
 			larcv::Image2D untagged = CollectedUntaggedPixels( img_v.at(p), thrumu_v.at(p), stopmu_v.at(p) );
 			untagged_v.emplace_back( std::move(untagged) );
 		}

 		// now we cluster the remaining hits
 		std::vector<untagged_cluster_info_t> untagged_cluster_info_v;
 		for (size_t p=0; p<untagged_v.size(); p++) {
 			untagged_cluster_info_t plane_cluster;
 			plane_cluster.pixels = dbscan::extractPointsFromImage( untagged_v.at(p), m_config.pixel_threshold );
 			dbscan::DBSCANAlgo algo;
 			plane_cluster.output = algo.scan( plane_cluster.pixels, 5, 10.0 );
 			untagged_cluster_info_v.emplace_back( std::move(plane_cluster) );
 		}

 		// with the clusters in hand (or in vector), we filter out the clusters we find interesting and locate extrema points.
 		// the extrema points serve as key points to use to try and match the cluster across planes

 	}

 	larcv::Image2D ContainedROI::CollectedUntaggedPixels( const larcv::Image2D& img, const larcv::Image2D& thrumu, const larcv::Image2D& stopmu) {

 		const larcv::ImageMeta& meta = img.meta();
 		larcv::Image2D untagged( meta );
 		untagged.paint(0.0);

 		for (size_t r=0; r<meta.rows(); r++) {
 			for (size_t c=0; c<meta.cols(); c++) {
 				if ( img.pixel(r,c) > m_config.pixel_threshold && (thrumu.pixel(r,c)==0 && stopmu.pixel(r,c)==0 ) ) {
 					untagged.set_pixel(r,c,1);
 				}
 			}
 		}

 		return untagged;
 	}

 	std::vector< ContainedROI::analyzed_cluster_t > ContainedROI::AnalyzeClusters( const untagged_cluster_info_t& clusters_info ) {

 		std::vector< ContainedROI::analyzed_cluster_t > interesting_clusters;

 		for (size_t idx_cluster=0; idx_cluster<clusters_info.output.clusters.size(); idx_cluster++ ) {

 			const dbscan::dbCluster& cluster = clusters_info.output.clusters.at(idx_cluster);

 			if ( (int)cluster.size() < m_config.min_cluster_size ) continue;

 			// we scan the hits, defining he extrema in row and column
 			std::vector< std::array<int,2> > extrema_pixels(4); // { highest in col, lowest in col, highest in row lowest in row }
 			std::vector< float > pixel_means(2,0.0); // mean in {col,row} dimensions
 			for (size_t ihit=0; ihit<cluster.size(); ihit++ ) {
 				int idx_hit = cluster.at(ihit);
 				const std::vector<double>& hit = clusters_info.pixels.at(idx_hit);

 				if ( ihit==0 ) {
 					// first hit, we fill the extrema store
 					for (size_t h=0; h<4; h++) {
 						extrema_pixels[h][0] = (int)hit[0];
	 					extrema_pixels[h][1] = (int)hit[1];
	 				}
 				}
 				else {
	 				// highest/lowest in col
	 				if ( extrema_pixels[0][0]<hit[0]) {
	 					extrema_pixels[0][0] = (int)hit[0];
	 					extrema_pixels[0][1] = (int)hit[1];
	 				}
	 				else if ( extrema_pixels[1][0]>hit[0] ) {
						extrema_pixels[1][0] = (int)hit[0];
	 					extrema_pixels[1][1] = (int)hit[1];
					}
	 				// highest/lowest in row
	 				if ( extrema_pixels[2][1]<hit[1]) {
	 					extrema_pixels[2][0] = (int)hit[0];
	 					extrema_pixels[2][1] = (int)hit[1];
	 				}
	 				else if ( extrema_pixels[3][1]>hit[1] ) {
						extrema_pixels[3][0] = (int)hit[0];
	 					extrema_pixels[3][1] = (int)hit[1];
					}
				}//end of ihit>0

				// find average of col and row
				for (size_t v=0; v<2; v++)
					pixel_means[v] += (float)hit[v];
 			}

 			// fill the output
 			for (size_t v=0; v<2; v++)
 				pixel_means[v] /= (float)cluster.size();

 			analyzed_cluster_t ana_cluster;
 			ana_cluster.cluster_idx = idx_cluster;
 			ana_cluster.mean = pixel_means;
 			for (size_t h=0; h<4; h++) {
				larcv::Pixel2D pix(extrema_pixels[h][0],extrema_pixels[h][1]);
				ana_cluster.extrema_pts.emplace_back( std::move(pix) );
 			}

 			interesting_clusters.emplace_back( std::move(ana_cluster) );
 		}//end of cluster loop

 		return interesting_clusters;

 	}


}