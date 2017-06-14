#include "ContainedROI.h"
#include <array>
#include <sstream>

// larcv
#include "UBWireTool/UBWireTool.h"

#ifndef __CINT__
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif
#endif

namespace larlitecv {


	ContainedROI::ContainedROI( const ContainedROIConfig& config ) {
		m_config = config;
		m_run = m_subrun = m_event = m_entry = -1;
	} 


 	std::vector<larcv::ROI> ContainedROI::SelectROIs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& thrumu_v,
    	const std::vector<larcv::Image2D>& stopmu_v, const std::vector<larcv::Image2D>& badch_v, 
    	std::vector< std::vector<larcv::Pixel2DCluster> >& matched_pixel_clusters ) {

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
 		m_untagged_cluster_info_v.clear();
 		for (size_t p=0; p<untagged_v.size(); p++) {
 			untagged_cluster_info_t plane_cluster;
 			plane_cluster.pixels = dbscan::extractPointsFromImage( untagged_v.at(p), 0.5 );
 			dbscan::DBSCANAlgo algo;
 			plane_cluster.output = algo.scan( plane_cluster.pixels, 5, 10.0 );
 			std::cout << "plane " << p << " has " << plane_cluster.output.clusters.size() << std::endl;
 			m_untagged_cluster_info_v.emplace_back( std::move(plane_cluster) );
 		}

#ifdef USE_OPENCV
 		std::vector< cv::Mat > cvimgs_v;
#endif

 		// with the clusters in hand (or in vector), we filter out the clusters we find interesting and locate extrema points.
 		// the extrema points serve as key points to use to try and match the cluster across planes
 		std::vector< std::vector<analyzed_cluster_t> > analyzed_clusters;
 		for (size_t p=0; p<untagged_v.size(); p++) {
			std::vector< ContainedROI::analyzed_cluster_t > plane_clusters = AnalyzeClusters( m_untagged_cluster_info_v.at(p), img_v.at(p) );
			std::cout << "plane " << p << " has " << plane_clusters.size() << " untagged clusters." << std::endl;
			const larcv::ImageMeta& meta = untagged_v.at(p).meta();
#ifdef USE_OPENCV
			cv::Mat imgmat = larcv::as_mat_greyscale2bgr( img_v.at(p), 5.0, 50.0 );
#endif

			for ( size_t icluster=0; icluster<plane_clusters.size(); icluster++ ) {
				const analyzed_cluster_t& cluster = plane_clusters.at(icluster);
				/*
				std::cout << " cluster #" << icluster << " index=" << cluster.cluster_idx
					<< " num hits=" << m_untagged_cluster_info_v.at(p).output.clusters.at(cluster.cluster_idx).size()
					<< " charge=" << cluster.total_charge
				  << " mean pixel=(" << cluster.mean[0] << "," << cluster.mean[1] << ") "
				  << " largest col=("  << cluster.extrema_pts.at(0).X() << "," << meta.pos_y(cluster.extrema_pts.at(0).Y()) << ") "
					<< " smallest col=(" << cluster.extrema_pts.at(1).X() << "," << meta.pos_y(cluster.extrema_pts.at(1).Y()) << ") "
				  << " largest row=("  << cluster.extrema_pts.at(2).X() << "," << meta.pos_y(cluster.extrema_pts.at(2).Y()) << ") "
				  << " smallest row=(" << cluster.extrema_pts.at(3).X() << "," << meta.pos_y(cluster.extrema_pts.at(3).Y()) << ") "
				  << std::endl;
				  */
#ifdef USE_OPENCV
				if (cluster.total_charge>m_config.min_cluster_plane_charge.at(p)) {
					cv::circle(imgmat,cv::Point(cluster.extrema_pts.at(0).X(),cluster.extrema_pts.at(0).Y()), 5, cv::Scalar(0,255,0),-1);
					cv::circle(imgmat,cv::Point(cluster.extrema_pts.at(1).X(),cluster.extrema_pts.at(1).Y()), 5, cv::Scalar(255,255,0),-1);
					cv::circle(imgmat,cv::Point(cluster.extrema_pts.at(2).X(),cluster.extrema_pts.at(2).Y()), 5, cv::Scalar(0,255,255),-1);
					cv::circle(imgmat,cv::Point(cluster.extrema_pts.at(3).X(),cluster.extrema_pts.at(3).Y()), 5, cv::Scalar(0,0,255),-1);
				}
#endif
			}

#ifdef USE_OPENCV
			cvimgs_v.emplace_back( std::move(imgmat) );
#endif

			analyzed_clusters.emplace_back( std::move(plane_clusters) );
		}

		// with list of clusters on each plane, we now make likelihood values for different combinations
		std::vector< std::vector<analyzed_cluster_t> > combinations;
		std::map<int,float> combo_likelihoods;
		std::map<int,float> combo_best_triarea;
		MatchClusters( analyzed_clusters.at(0), analyzed_clusters.at(1), analyzed_clusters.at(2), img_v, combinations, combo_likelihoods, combo_best_triarea );
		std::cout << "Number of cluster combinations: " << combinations.size() << std::endl;
		std::vector< std::pair<int,float> > ordered_combos;
		for ( size_t icombo=0; icombo<combinations.size(); icombo++ ) {
			const std::vector< analyzed_cluster_t >& combo = combinations.at(icombo);
			float ll = combo_likelihoods[icombo];
			//std::cout << "Combo #" << icombo << " (" << ")" << " ll=" << ll << std::endl;
			ordered_combos.push_back( std::pair<int,float>(icombo,ll) );
		}

		// sort the combinations
		struct combo_sorter {
			bool operator()( std::pair<int,float> lhs, std::pair<int,float> rhs ) {
				if ( lhs.second<rhs.second) return true;
				return false;
			};
		} my_sorter;
		std::sort( ordered_combos.begin(), ordered_combos.end(), my_sorter );

		// evaluate and save ROIs
		std::vector<larcv::ROI> output;

		// clear pixel clusters
		matched_pixel_clusters.clear();

		int ncombos = 0;
		std::vector< std::set<int> > used_clusters;
		for ( auto const& combo_idx : ordered_combos ) {
			const std::vector<analyzed_cluster_t>& combo = combinations.at(combo_idx.first);
			

			// acceptance conditions...			
			// we know this will be in adequate, but we use some threshold for the likelihood
			if ( combo_idx.second < m_config.roi_likelihood_cutoff ) {
				// we accept, make an roi
				larcv::ROI untagged_roi;
				std::vector<larcv::Pixel2DCluster> matched_cluster;

				for ( size_t p=0; p<3; p++ ) {
					const larcv::ImageMeta& meta = img_v.at(p).meta();
					const analyzed_cluster_t& cluster = combo.at(p);
					// use extrema to define bounding box
					int drow = cluster.extrema_pts.at(2).Y() - cluster.extrema_pts.at(3).Y();
					int dcol = cluster.extrema_pts.at(0).Y() - cluster.extrema_pts.at(1).Y();
					float min_x = meta.pos_x( cluster.extrema_pts.at(1).X() );
					float max_x = meta.pos_x( cluster.extrema_pts.at(0).X() );
					float min_y = meta.pos_y( cluster.extrema_pts.at(2).Y() );
					float max_y = meta.pos_y( cluster.extrema_pts.at(3).Y() );
					larcv::ImageMeta bb( fabs(max_x-min_x), fabs(max_y-min_y), drow, dcol, min_x, max_y, (larcv::PlaneID_t)p );
					untagged_roi.AppendBB( bb );

					larcv::Pixel2DCluster pix_cluster;
					const dbscan::dbCluster& hitcluster = m_untagged_cluster_info_v.at(p).output.clusters.at( cluster.cluster_idx );
					for ( auto &hitidx : hitcluster ) {
						larcv::Pixel2D pixhit( m_untagged_cluster_info_v.at(p).pixels.at(hitidx)[0], m_untagged_cluster_info_v.at(p).pixels.at(hitidx)[1] );
						pix_cluster += pixhit;
					}
					matched_cluster.emplace_back( std::move(pix_cluster) );
				}

				/*
				std::cout << "Combo #" << ncombos 
					<< "(" << combo.at(0).cluster_idx << "," << combo.at(1).cluster_idx << "," << combo.at(2).cluster_idx << ") "
					<< " ll=" << combo_idx.second 
					<< " ll(triarea)=" << combo_best_triarea[combo_idx.first]
					<< " y-plane bounds=[" << untagged_roi.BB(2).min_x() << "-" << untagged_roi.BB(2).max_x() << "]"
					<< std::endl;
					*/

				output.emplace_back( std::move(untagged_roi) );
				matched_pixel_clusters.emplace_back( std::move(matched_cluster) );
			}

			ncombos++;
			if ( ncombos>=m_config.max_number_rois )
				break;
			if ( ncombos>m_config.max_number_rois && combo_idx.second > m_config.max_roi_likelihood )
				break;
		}

#ifdef USE_OPENCV
		for ( size_t p=0; p<3; p++ ) {
			std::stringstream ss;
			ss << "contained_roi_";
			if ( m_run!=-1 &&m_subrun!=-1 && m_event!=-1 && m_entry!=-1 ) {
				ss << m_entry << "_" << m_run << "_" << m_subrun << "_" << m_event << "_";
			}
			ss << "p" << p << ".jpg";

			// add roi's
			for ( auto &roi : output ) {
				larcv::draw_bb( cvimgs_v.at(p), img_v.at(p).meta(), roi.BB().at(p), 0, 200, 0, 2 );
			}

			if ( m_config.draw_truth_roi && m_truth_roi.BB().size()>0 ) {
				larcv::draw_bb( cvimgs_v.at(p), img_v.at(p).meta(), m_truth_roi.BB().at(p), 200, 0, 0, 2 );
			}
			//cv::imwrite( ss.str(), cvimgs_v.at(p) );
		}
#endif

		return output;

 	}

 	larcv::Image2D ContainedROI::CollectedUntaggedPixels( const larcv::Image2D& img, const larcv::Image2D& thrumu, const larcv::Image2D& stopmu) {

 		const larcv::ImageMeta& meta = img.meta();
 		larcv::Image2D untagged( meta );
 		untagged.paint(0.0);

 		for (size_t r=0; r<meta.rows(); r++) {
 			for (size_t c=0; c<meta.cols(); c++) {
 				if ( img.pixel(r,c) > m_config.pixel_threshold && (thrumu.pixel(r,c)<0.5 && stopmu.pixel(r,c)<0.5 ) ) {
 					untagged.set_pixel(r,c,1);
 				}
 			}
 		}

 		return untagged;
 	}

 	std::vector< ContainedROI::analyzed_cluster_t > ContainedROI::AnalyzeClusters( const untagged_cluster_info_t& clusters_info, const larcv::Image2D& img ) {

 		std::vector< ContainedROI::analyzed_cluster_t > interesting_clusters;

 		for (size_t idx_cluster=0; idx_cluster<clusters_info.output.clusters.size(); idx_cluster++ ) {

 			const dbscan::dbCluster& cluster = clusters_info.output.clusters.at(idx_cluster);

 			if ( (int)cluster.size() < m_config.min_cluster_size ) continue;

 			// we scan the hits, defining he extrema in row and column
 			std::vector< std::array<int,2> > extrema_pixels(4); // { highest in col, lowest in col, highest in row lowest in row }
 			std::vector< float > pixel_means(2,0.0); // mean in {col,row} dimensions
 			float total_charge = 0.;
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

				// total the charge
				total_charge += img.pixel( (int)hit[1], (int)hit[0] );
 			}

 			// fill the output
 			for (size_t v=0; v<2; v++)
 				pixel_means[v] /= (float)cluster.size();

 			analyzed_cluster_t ana_cluster;
 			ana_cluster.cluster_idx = idx_cluster;
 			ana_cluster.mean = pixel_means;
 			ana_cluster.total_charge = total_charge;
 			for (size_t h=0; h<4; h++) {
				larcv::Pixel2D pix(extrema_pixels[h][0],extrema_pixels[h][1]);
				ana_cluster.extrema_pts.emplace_back( std::move(pix) );
 			}

 			interesting_clusters.emplace_back( std::move(ana_cluster) );
 		}//end of cluster loop

 		return interesting_clusters;

 	}

 	float ContainedROI::CalculateMatchLikelihood( const std::vector<analyzed_cluster_t>& clusters, const std::vector<larcv::Image2D>& img_v, 
 		std::vector<float>& ll_components ) {
 		// we build a matching score based on the following quantities
 		// 1) charge
 		// 2) timing of start and end of cluster
 		// 3) 3D consistency 

 		int nplanes = (int)clusters.size();

 		typedef std::pair<int,int> plane_combo_t;

 		// 1) charge and end timings
 		std::map< plane_combo_t, float> charge_diffs;
 		std::map< plane_combo_t, float> start_diffs;
	 	std::map< plane_combo_t, float> end_diffs;

 		for (int a=0; a<nplanes; a++) {
 			for (int b=a+1; b<nplanes; b++) {

 				const analyzed_cluster_t& cluster_a = clusters.at(a);
 				const analyzed_cluster_t& cluster_b = clusters.at(b);
 				const larcv::ImageMeta& meta_a = img_v.at(a).meta();
 				const larcv::ImageMeta& meta_b = img_v.at(b).meta();

 				// charge
		 		float charge_diff = cluster_a.total_charge - cluster_b.total_charge;
		 		float ave_charge = 0.5*(cluster_a.total_charge + cluster_b.total_charge);
		 		charge_diff /= ave_charge; // needs to be with respect to total charge
			  plane_combo_t q_combo(a,b);
		 		charge_diffs.insert( std::make_pair< plane_combo_t, float >( std::move(q_combo), std::move(charge_diff) ) );

		 		// end timings
		 		plane_combo_t start_combo(a,b);
		 		float start_diff = meta_a.pos_y(cluster_a.extrema_pts.at(3).Y()) - meta_b.pos_y(cluster_b.extrema_pts.at(3).Y());
		 		start_diffs.insert( std::make_pair< plane_combo_t, float >(std::move(start_combo),std::move(start_diff) ) );

		 		plane_combo_t end_combo(a,b);
 				float end_diff   = meta_a.pos_y(cluster_a.extrema_pts.at(2).Y()) - meta_b.pos_y(cluster_b.extrema_pts.at(2).Y());
 				end_diffs.insert(   std::make_pair< plane_combo_t, float >(std::move(end_combo),std::move(end_diff)) );

 			}
 		}

 		// 3) 3D consistency of extrema 
 		std::vector<float> extrema_triarea(4,1.0e6);
 		for (int ipt=0; ipt<4; ipt++) {
 			// get wid of all three planes
 			std::vector<int> wid(3,0);
 			for (int p=0; p<nplanes; p++) {
	 			wid[p] = img_v.at(p).meta().pos_x( clusters.at(p).extrema_pts.at(ipt).X() );
	 		}
	 		double triarea = 0;
	 		int crosses = 0;
	 		std::vector<float> intersection_yz;
	 		larcv::UBWireTool::wireIntersection( wid, intersection_yz, triarea, crosses );
	 		if ( crosses==1 ) {
	 			extrema_triarea[ipt] = triarea;
	 		}
 		}

 		// calculate the likelihood 
 		float likelihood = 0.;
 		for ( auto& it_charge_diff : charge_diffs ) {
 			float cdiff = it_charge_diff.second/m_config.charge_diff_sigma;
 			float ll_charge = m_config.charge_diff_weight*0.5*cdiff*cdiff;
 			likelihood += ll_charge;
 			ll_components.push_back( ll_charge );
 		}

 		for ( auto& it_start_diff : start_diffs ) {
 			float sdiff = it_start_diff.second/m_config.time_boundary_diff_sigma;
 			float ll_start = m_config.time_boundary_diff_weight*0.5*sdiff*sdiff;
 			likelihood += ll_start;
 			ll_components.push_back( ll_start );
 		}

 		for ( auto& it_end_diff : end_diffs ) {
 			float ediff = it_end_diff.second/m_config.time_boundary_diff_sigma;
 			float ll_end = m_config.time_boundary_diff_weight*0.5*ediff*ediff; 
 			likelihood += ll_end;
 			ll_components.push_back( ll_end );
 		}

 		float triarea_likelihood = 1.0e6;
 		float best_triarea = 1.0e6;
 		for ( auto& triarea : extrema_triarea ) {
 			if ( triarea >= 1.0e6 ) continue;
 			float tdiff = triarea/m_config.triarea_sigma;
 			float ll_triarea = m_config.triarea_weight*0.5*tdiff*tdiff;
 			if ( best_triarea > triarea ) {
	 			triarea_likelihood = ll_triarea;
 				best_triarea = triarea;
 			}	
 		}
 		likelihood += triarea_likelihood;
 		ll_components.push_back(best_triarea);

 		return likelihood;
 	}


	void ContainedROI::MatchClusters( const std::vector<analyzed_cluster_t>& uplane, const std::vector<analyzed_cluster_t>& vplane, const std::vector<analyzed_cluster_t>& yplane, 
			const std::vector<larcv::Image2D>& img_v,  std::vector< std::vector<analyzed_cluster_t> >& cluster_combos, 
			std::map< int, float >& combo_likelihoods, std::map<int,float>& combo_best_triarea ) {

		cluster_combos.clear();
		combo_likelihoods.clear();

		for ( auto const& ucluster : uplane ) { 
			if ( ucluster.total_charge < m_config.min_cluster_plane_charge.at(0) ) continue;
			for ( auto const& vcluster: vplane ) {
				if ( vcluster.total_charge < m_config.min_cluster_plane_charge.at(1) ) continue;
				for ( auto const& ycluster : yplane ) {
					if ( ycluster.total_charge < m_config.min_cluster_plane_charge.at(2) ) continue;

					// score this combination
					std::vector< analyzed_cluster_t > candidate_combo;
					candidate_combo.push_back( ucluster );
					candidate_combo.push_back( vcluster );
					candidate_combo.push_back( ycluster );
					std::vector<float> ll_components;
					float ll = CalculateMatchLikelihood( candidate_combo, img_v, ll_components );
					float ll_best_tri = ll_components.back();
					cluster_combos.push_back( std::move(candidate_combo) );
					int combo_idx = (int)cluster_combos.size() - 1;
					combo_likelihoods.insert( std::make_pair<int,float>( std::move(combo_idx), std::move(ll) ) );
					combo_best_triarea.insert( std::make_pair<int,float>( std::move(combo_idx), std::move(ll_best_tri) ) );
				}
			}
		}
	}

}
