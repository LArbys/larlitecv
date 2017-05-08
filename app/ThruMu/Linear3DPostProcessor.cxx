#include "Linear3DPostProcessor.h"

// larlite/UserDev/BasicTool/GeoAlgo
#include "BasicTool/GeoAlgo/GeoAlgo.h"
#include "BasicTool/GeoAlgo/GeoLineSegment.h"

//larlitecv/app
#include "TaggerTypes/BoundarySpacePoint.h"

namespace larlitecv {

	std::vector< BMTrackCluster3D > Linear3DPostProcessor::process( std::vector< BMTrackCluster3D >& tracks_v ) {

		geoalgo::GeoAlgo algo;

		// we pair-wise check for subsets
		const size_t ntracks = tracks_v.size();
		std::vector<bool> track_excluded( ntracks, false );

		for (size_t itrack_a=0; itrack_a<ntracks; itrack_a++) {
			if ( track_excluded.at(itrack_a) ) continue;

			BMTrackCluster3D& track_a = tracks_v.at(itrack_a);
			geoalgo::LineSegment seg_a( track_a.path3d.front()[0], track_a.path3d.front()[1], track_a.path3d.front()[2],
																	track_a.path3d.back()[0],  track_a.path3d.back()[1],  track_a.path3d.back()[2] );

			for (size_t itrack_b=itrack_a+1; itrack_b<ntracks; itrack_b++ ) {
				if ( track_excluded.at(itrack_b) ) continue;

				BMTrackCluster3D& track_b = tracks_v.at(itrack_b);
				geoalgo::LineSegment seg_b( track_b.path3d.front()[0], track_b.path3d.front()[1], track_b.path3d.front()[2],
																		track_b.path3d.back()[0],  track_b.path3d.back()[1],  track_b.path3d.back()[2] );

				// same start or end point
				BoundarySpacePoint* sppt[2][2]    = { { &track_a.start_endpts, &track_a.end_endpts},
																							{ &track_b.start_endpts, &track_b.end_endpts } };
				std::vector<double>* endpts[2][2] = { { &track_a.path3d.front(), &track_a.path3d.back() },
																							{ &track_b.path3d.front(), &track_b.path3d.back() } };
				bool issame[2] = { false, false };
				for (int ipt=0; ipt<2; ipt++) {
					float dist = 0.;
					for (int i=0; i<3; i++)	{
						dist += (endpts[ipt][1]->at(i) - endpts[ipt][0]->at(i))*(endpts[ipt][1]->at(i) - endpts[ipt][0]->at(i));
					}
					dist = sqrt(dist);
					if ( dist<1.0e-5 && sppt[ipt][0]->type()==sppt[ipt][1]->type() ) {
						issame[ipt] = true;
					}
				}

				float pair_cosine = fabs( seg_a.Dir().Dot( seg_b.Dir() )/( seg_a.Dir().Length()*seg_b.Dir().Length() ) );

				std::cout << "Criteria #1: Comparing track[" << itrack_a  << "] and track[" << itrack_b << "]" << std::endl;
				std::cout << "  are same: start=" << issame[0] << " end=" << issame[1] << std::endl;
				std::cout << "  cosine: " << pair_cosine << std::endl;

				// critera #1: if one end is the same and in the same direction
				if ( (issame[0] || issame[1]) && pair_cosine>0.8 ) {
					// we exclude the shortest one
					if ( seg_a.Dir().Length()<seg_b.Dir().Length() ) {
						track_excluded.at(itrack_a) = true;
						break; // track-a is excluded. no need to check it against others
					}
					else {
						track_excluded.at(itrack_b) = true;
					}
					continue;
				}

				// neither end point is the same

				// if end points of shortest track is close to longest track, it's a duplicate
				geoalgo::LineSegment* shorter_lineseg = &seg_b;
				geoalgo::LineSegment* longer_lineseg  = &seg_a;
				size_t shorter_seg_idx = itrack_b;
				//size_t longer_seg_idx  = itrack_a;
				if ( seg_a.Dir().Length()<seg_b.Dir().Length() ) {
					// flip it
					shorter_lineseg = &seg_a;
					longer_lineseg  = &seg_b;
					shorter_seg_idx = itrack_a;
					//longer_seg_idx = itrack_b;
				}
				geoalgo::Line longer_line( longer_lineseg->Start(), longer_lineseg->End() );
				float start_dist = algo.SqDist( shorter_lineseg->Start(), longer_line );
				float end_dist   = algo.SqDist( shorter_lineseg->End(),   longer_line );

				std::cout << "Criteria #2: Comparing track[" << itrack_a  << "] and track[" << itrack_b << "] "
									<< "  have start and end close to tracks. " << std::endl;
				std::cout << "  startdist=" << start_dist << std::endl;
				std::cout << "  enddist=" << end_dist << std::endl;

				// criterion #2
				if ( start_dist<20.0 && end_dist<20.0 ) {

					// criterion #2A: shortest track is contained inside longer one


					// criterion #2B: does short track end extend past longer one, if so, we extend the longer one

					// ok, we're done with one of these tracks
					track_excluded.at( shorter_seg_idx ) = true;
					if ( shorter_seg_idx==itrack_a ) 
						break;
					else
						continue;
				}

			}
		}

		// collect the output
		std::vector< BMTrackCluster3D > output_tracks;

		for ( size_t itrack=0; itrack<track_excluded.size(); itrack++ )	{
			if ( track_excluded.at(itrack) ) 
				continue;

			output_tracks.emplace_back( std::move(tracks_v.at(itrack)) );
		}

		tracks_v.clear();

		return output_tracks;

	}

}
