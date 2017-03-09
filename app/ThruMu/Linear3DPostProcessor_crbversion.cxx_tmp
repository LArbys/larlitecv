#include "Linear3DPostProcessor.h"

// larlite/UserDev/BasicTool/GeoAlgo                                                                                                              
#include "GeoAlgo.h"
#include "GeoLineSegment.h"

//larlitecv/app/ThruMu                                                                                                                         
#include "BoundarySpacePoint.h"

// Include the '<vector>' class for all of its functionality.
#include <vector>


namespace larlitecv {

  // Declare the function that will take care of processing the linear tracks after the fact.
  // Call the function the 'process' function.
  // Input: (1) 'tracks_v' - This is the collection of tracks that were initially tagged by the '3DLinearFitter' and the 'AStar3D' Algo.
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
	  for (int i=0; i<3; i++){
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
	size_t longer_seg_idx  = itrack_a;
	if ( seg_a.Dir().Length()<seg_b.Dir().Length() ) {
	  // flip it
	  shorter_lineseg = &seg_a;
	  longer_lineseg  = &seg_b;
	  shorter_seg_idx = itrack_a;
	  longer_seg_idx = itrack_b;
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
	  // This can be accomplished by a twofold process
	  // (1) First, ensure that that tracks from which both segments are derived are moving in the same direction (use the check on 'cos' used above).
	  // (2) Second, make sure that every point on 'shorter_lineseg' is a minimum distance away from a point on 'longer_lineseg'.
	    
	  // Ensure that the tracks are moving in the same direction with an outer loop. 
	  if (pair_cosine > 0.8) {

	    // Use the 'pointProximity' function to find if every point on the shorter line segment is a maximum distance of 5 cm away from the longer line segment.
	    // If this is the case, then you can eliminate the shorter track because it is completely enclosed by the longer track.
	    if (pointProximity(longer_line, shorter_lineseg, 10, 5.0) == true) {
	            
	      track_excluded.at( shorter_seg_idx ) = true;

	      // You can break from the inner loop, because you have already removed the inner track.
	      if ( shorter_seg_idx==itrack_a )
		break;
	            
	      // Otherwise, continue on with the next iteration of the loop.
	      else
		continue;

	    }

	    // Use the 'longerTrackExtensionCoordinates' to find out how you should extend the track in question.
	    // I can extend the track in question now.
	    // The output will be a vector of the coordinates of the two new endpoints of the new track. (a 2 * 3 vector)
	    std::vector < std::vector <double> > endpoint_coordinates_of_new_track = longerTrackExtensionCoordinates(longer_line.Pt1(), longer_line.Pt2(), shorter_lineseg->Start(), shorter_lineseg->End());

	    // If the two tracks do not extend past one another, then the every entry of the 'endpoint_coordinates_of_new_track' vector will be equal to 0.
	    // Check the first entry of the first set of coordinates to see if it is 0, and if so continue on to the next point in the loop.  These two tracks satisfy none of the conditions for eliminating one of them (they don't have the same endpoints, and the endpoints are close by but not every point on the line segment is close to every point on the track and the line segment does not extend out from the longer line.
	    if (endpoint_coordinates_of_new_track[0][0] > -0.001 && endpoint_coordinates_of_new_track[0][0] < 0.001)
	      continue;
	    
	    // Consider the case now that the coordinate values within 'endpoint_coordinates_of_new_track' are nonzero.
	    else if (endpoint_coordinates_of_new_track[0][0] < -0.001 || endpoint_coordinates_of_new_track[0][0] > 0.001) {
	      
	      // Adjust 'tracks_v' by using the 'concatenatingTwoTrackParts' function.
	      concatenatingTwoTrackParts(tracks_v, itrack_a, itrack_b, track_a, track_b, endpoint_coordinates_of_new_track);

	          
	      // Once you have made it this far, you can exclude the smaller track segment from the analysis.

	      track_excluded.at( shorter_seg_idx ) = true;

	      // You can break from the inner loop, because you have already removed the inner track.
	      if ( shorter_seg_idx==itrack_a )
		break;
	            
	      // Otherwise, continue on with the next iteration of the loop.
	      else
		continue;
	            
	    } // End of the 'else if' considering that one track extends past the other

	  } // End of the case in which the two tracks almost point in the same direction

	} // End of the case where the two track endpoints are close

      } // End of the inner loop over the tracks

    } // End of the outer loop over the tracks

    // collect the output
    std::vector< BMTrackCluster3D > output_tracks;

    for ( size_t itrack=0; itrack<track_excluded.size(); itrack++ ){
      if ( track_excluded.at(itrack) ) 
	continue;

      output_tracks.emplace_back( std::move(tracks_v.at(itrack)) );
    }

    tracks_v.clear();

    return output_tracks;

  } // End of the 'process' function


  // A function that finds if the points on the line segment are close to those on the line.                                                                                
  // Input values: (1) 'long_line' - an object of type 'Line' that contains the initial point and the final point on the longer of the two track-like entities.                
  //               (2) 'short_linesegment' - a line segment that is the shorter of the two track-like entities.                                                               
  //               (3) 'resolution' - The number of line segments that you want on both 'short_linesegment' and 'long_line' when comparing to see how close they are to one another.
  //               (4) 'maximum_dist' - The maximum distance that a point on 'short_linesegment' can be from any point on 'long_line' while still considering the two overlapping.
  bool Linear3DPostProcessor::pointProximity(geoalgo::Line long_line, geoalgo::LineSegment* short_linesegment, int resolution, double maximum_dist) {

    // Unpack the 'long_line' endpoint coordinates into a vector.                                                                                                             
    const std::vector < double > long_line_first_endpoint_coords(long_line.Pt1());            
    
    const std::vector < double > long_line_second_endpoint_coords(long_line.Pt2());

    // Unpack the 'short_linesegment' endpoints in a vector.
    const std::vector < double > shorter_linesegment_first_endpoint_coords(short_linesegment->Start());
    
    const std::vector < double > shorter_linesegment_second_endpoint_coords(short_linesegment->End()); 

    // Find the dimension with the greatest extent for both the line and the short line segment.                                                                            
    // The index notation that I will use for the dimensions with respect to one another is '0' - x, '1' - y, '2' - z.                                                    
    // line.                                                                                                                                                                   
    int    line_longest_dim_idx              = longestDimFinder(long_line_first_endpoint_coords, long_line_second_endpoint_coords);
    int    short_linesegment_longest_dim_idx = longestDimFinder(shorter_linesegment_second_endpoint_coords, shorter_linesegment_second_endpoint_coords);

    // Declare a vector for the two indices that do NOT correspond to the dimension with the furthest spread in its dimension.                                              
    std::vector <int> line_shorter_dim_idx_vector;
    std::vector <int> short_linesegment_shorter_dim_idx_vector;

    // Now, you can use a loop to find out which dimension of proximal points are smaller than the largest.                 
    for (size_t dim_iter = 0; dim_iter < 3; dim_iter++) {

      // Check to see if this dimension is not equal to 'line_longest_dim_idx'.                                                              
      if (dim_iter != line_longest_dim_idx)
        // If the conditional is satisfied, then append 'dim_iter' onto 'line_shorter_dim_idx_vector'.                                                                  
	line_shorter_dim_idx_vector.push_back(dim_iter);

      // Check to see if this dimension is not equal to 'short_linesegment_longest_dim_idx'.                                                                                  
      if (dim_iter != short_linesegment_longest_dim_idx)
        // If the conditional is satisfied, then append 'dim_iter' onto 'short_linesegment_shorter_dim_idx_vector'.                                                              
	short_linesegment_shorter_dim_idx_vector.push_back(dim_iter);

    }

    // Generate the points on each of the lines using the points on the two track-like entities, the resolution, the '_shorter_dim_idx_vector' for each track-like entity, and the '_longest_dim_idx' value for each track-like entity.                                                                                                                    
    std::vector < std::vector<double> > long_line_points           = getTrackPoints(long_line_first_endpoint_coords, long_line_second_endpoint_coords, resolution, line_shorter_dim_idx_vector, line_longest_dim_idx);
    std::vector < std::vector<double> > short_linesegment_points   = getTrackPoints(shorter_linesegment_first_endpoint_coords, shorter_linesegment_second_endpoint_coords, resolution, short_linesegment_shorter_dim_idx_vector, short_linesegment_longest_dim_idx);

    // Declare a variable for if this point on the line segment is close to any of the points on the line.  Initialize it to 'True'.                         
    bool is_linesegment_point_close_to_point_on_line = false;

    // Start a loop over the points on the short linesegment as the outer loop.  Every point on the short linesegment will be compared.
    for (size_t segment_pt_iter = 0; segment_pt_iter < short_linesegment_points.size(); segment_pt_iter++) {

      // Reset 'is_linesegment_point_close_to_point_on_line' to 'False'.                                                                                                       
      is_linesegment_point_close_to_point_on_line  = false;

      for (size_t line_pt_iter = 0; line_pt_iter < long_line_points.size(); line_pt_iter++) {

        if (sqrt((long_line_points[line_pt_iter][0] - short_linesegment_points[segment_pt_iter][0]) * (long_line_points[line_pt_iter][0] - short_linesegment_points[segment_pt_iter][0]) + (long_line_points[line_pt_iter][1] - short_linesegment_points[segment_pt_iter][1]) * (long_line_points[line_pt_iter][1] - short_linesegment_points[segment_pt_iter][1]) + (long_line_points[line_pt_iter][2] - short_linesegment_points[segment_pt_iter][2]) * (long_line_points[line_pt_iter][2] - short_linesegment_points[segment_pt_iter][2])) < maximum_dist) {

          // Set 'is_linesegment_point_close_to_point_on_line' to 'True'.                                                                                                   
	  is_linesegment_point_close_to_point_on_line  = true;

          // Break the loop.  It needn't continue any more for its purposes.                                                                                   
	  break;

	}
	
      }

      // Check to see the value of 'is_linesegment_point_close_to_point_on_line'.  If it is 'False', then return 'False' for the function.                           
      if (is_linesegment_point_close_to_point_on_line == false)
	return false;

}

    // If the sequence of 'for' loops above did not return 'false', then return 'true'.  Every point on the segment is close to every point on the line.           
    return true;

}


  // Use a function to determine if you should extend the longer linesegment.  If so, which points should mark the ends of the new track within the 'BMTTrackCluster3D' object in the 'track_v' vector.
  // Inputs: (1) 'longer_line_first_point_coords'          - The coordinates of the first point of the longer line of the two track-like objects.
  //         (2) 'longer_line_second_point_coords'         - The coordinates of the second point of the longer line of the two track-like objects.
  //         (3) 'shorter_linesegment_first_point_coords'  - The coordinates of the first point of the shorter line segment of the two track-like objects.
  //         (4) 'shorter_linesegment_second_point_coords' - The coordinates of the second point of the shorter line segment of the two track-like objects
  std::vector <std::vector <double> > longerTrackExtensionCoordinates(std::vector<double> longer_line_first_point_coords, std::vector<double> longer_line_second_point_coords, std::vector<double> shorter_linesegment_first_point_coords, std::vector<double> shorter_linesegment_second_point_coords) {

    // I use the assumption that the track extends in a straight line for this problem, meaning that each of the track's coordinates are continuously increasing or continuously decreasing for their entire path.                                                                                                                       

    // Find out if the line segment's coordinates are increasing along its path.                                                                     
    // Initialize variables for this first.                                                                                                          
    bool x_is_increasing_along_path_shorter_linesegment = true;
    bool y_is_increasing_along_path_shorter_linesegment = true;
    bool z_is_increasing_along_path_shorter_linesegment = true;

    // Reset each of these variables if they were not properly initialized.                                                                          
    // 'x_is_increasing_along_path_shorter_linesegment'                                                                                              
    if (shorter_linesegment_first_point_coords[0] > shorter_linesegment_second_point_coords[0])
      x_is_increasing_along_path_shorter_linesegment = false;

    // 'y_is_increasing_along_path_shorter_linesegment'                                                                                              
    if (shorter_linesegment_first_point_coords[1] > shorter_linesegment_second_point_coords[1])
      y_is_increasing_along_path_shorter_linesegment = false;

    // 'z_is_increasing_along_path_shorter_linesegment'                                                                                              
    if (shorter_linesegment_first_point_coords[2] > shorter_linesegment_second_point_coords[2])
      z_is_increasing_along_path_shorter_linesegment = false;


    // Find out if the longer line's coordinates are increasing along their path.                                                                    
    // Initialize variables for this first.                                                                                                          
    bool x_is_increasing_along_path_longer_line = true;
    bool y_is_increasing_along_path_longer_line = true;
    bool z_is_increasing_along_path_longer_line = true;

    // Reset each of these variables if they were not properly initialized.                                                                          
    // 'x_is_increasing_along_path_longer_line'                                                                                                      
    if (longer_line_first_point_coords[0] > longer_line_second_point_coords[0])
      x_is_increasing_along_path_longer_line = false;

    // 'y_is_increasing_along_path_longer_line'                                                                                                      
    if (longer_line_first_point_coords[1] > longer_line_second_point_coords[1])
      y_is_increasing_along_path_longer_line = false;

    // 'z_is_increasing_along_path_longer_line'                                                                                                      
    if (longer_line_first_point_coords[2] > longer_line_second_point_coords[2])
      z_is_increasing_along_path_longer_line = false;


    // I will solve this with a simple trick: I will make sure that the points that I am comparing between all sets of coordinates are the same and that the extremity requirement (that the track be greater than or less than a specific value) is satisfied.  Otherwise, I will return a vector that entirely consists of zeros, which is a sign that the extremity requirements were not satisfied.
    
    // Declare vectors of type 'int' for each of the variables.  The index system for these points is the following:
    // '0' if it is the first point on the track's trajectory that should be compared.
    // '1' if it is the second point on the track's trajectory that should be compared.
    // I will initialize these vectors to have length '2' and to consist of two 'int' values.

    // These vectors must all have the same value at the end, or else the shorter track is not an extension of the longer track.
    
    std::vector <int> x_points_to_compare(2, 0);
    std::vector <int> y_points_to_compare(2, 0);
    std::vector <int> z_points_to_compare(2, 0);

    // In a very clunky way, declare the vector that will be returned if none of the conditions are satisfied.
    std::vector <std::vector < double > > null_vector;
    null_vector.push_back(std::vector <double> (2, 0.0));
    null_vector.push_back(std::vector <double> (2, 0.0));

    // I will declare a 2D vector that I will call 'points_to_compare'.  This will take in '0' if I want to the compare the first point of the track, and it will take in '1' if I want to compare the second point of the track.

    // I will take the possible combination of 'true' and 'false' statements (all of them), and give to a function the variables that will have to be compared.

    // Use a series of 'if' loops to output the correct coordinates that should be compared.  Return the default if the condition is not satisfied.
    if (x_is_increasing_along_path_shorter_linesegment == true && x_is_increasing_along_path_longer_line == true) {

      
      // Compare the two starting x-coordinates to see if the 'shorter_linesegment' x-coordinate is less than that of the 'longer_line' x-coordinate.  Then the two x-coordinates to compare are the two starting coordinates, or [0, 0] in our notation.
      if (shorter_linesegment_first_point_coords[0] < longer_line_first_point_coords[0]) {

	// Set 'x_points_to_compare' to [0, 0]
	x_points_to_compare[0] = 0;
	x_points_to_compare[0] = 0;

      }

      // Compare the two ending x-coordinates to see if the 'longer_linesegment' x-coordinate is greater than that of the 'longer_line' x-coordinate.  Then the two x-coordinates to compare are the two ending coordinates, or [1, 1] in our notation.
      if (shorter_linesegment_second_point_coords[0] > longer_line_second_point_coords[0]) {

	// Set 'x_points_to_compare' to [1, 1]
	x_points_to_compare[0] = 1;
	x_points_to_compare[1] = 1;

      }

      // If neither of these conditions are satisfied, then the shorter linesegment does not extend past the longer track.  Just return a list of zeros, which does not compromise the type of the function and indicates that the shorter track does not extend past the longer track.
      else {

	return null_vector;

      }

    }

    // I will mute the comments now because of the great many similar loops.
    
    if (x_is_increasing_along_path_shorter_linesegment == true && x_is_increasing_along_path_longer_line == false) {
      
      if (shorter_linesegment_first_point_coords[0] < longer_line_second_point_coords[0]) {

	x_points_to_compare[0] = 0;
	x_points_to_compare[1] = 1;

      }

      if (shorter_linesegment_second_point_coords[0] > longer_line_first_point_coords[0]) {

	x_points_to_compare[0] = 1;
	x_points_to_compare[1] = 0;

      }

      else {

	return null_vector;

      }

    }

    if (x_is_increasing_along_path_shorter_linesegment == false && x_is_increasing_along_path_longer_line == true) {

      if (shorter_linesegment_second_point_coords[0] < longer_line_first_point_coords[0]) {

	x_points_to_compare[0] = 1;
	x_points_to_compare[1] = 0;

      }

      if (shorter_linesegment_first_point_coords[0] > longer_line_second_point_coords[0]) {

	x_points_to_compare[0] = 0;
	x_points_to_compare[1] = 1;

      }

      else {

	return null_vector;

      }

    }

    if (x_is_increasing_along_path_shorter_linesegment == false && x_is_increasing_along_path_longer_line == false) {

      if (shorter_linesegment_second_point_coords[0] < longer_line_second_point_coords[0]) {

	x_points_to_compare[0] = 1;
	x_points_to_compare[0] = 1;

      }

      if (shorter_linesegment_first_point_coords[0] > longer_line_first_point_coords[0]) {

	x_points_to_compare[0] = 0;
	x_points_to_compare[0] = 0;

      }
   
      else {

        return null_vector;

      }

    }

    // y-coordinates

    if (y_is_increasing_along_path_shorter_linesegment == true && y_is_increasing_along_path_longer_line == true) {

      if (shorter_linesegment_first_point_coords[1] < longer_line_first_point_coords[1]) {

	y_points_to_compare[0] = 0;
	y_points_to_compare[1] = 0;

      }

      if (shorter_linesegment_second_point_coords[1] > longer_line_second_point_coords[1]) {

	y_points_to_compare[0] = 1;
	y_points_to_compare[1] = 1;

      }

      else {

	return null_vector;

      }

    }


    if (y_is_increasing_along_path_shorter_linesegment == true && y_is_increasing_along_path_longer_line == false) {

      if (shorter_linesegment_first_point_coords[1] < longer_line_second_point_coords[1]) {

	y_points_to_compare[0] = 0;
	y_points_to_compare[1] = 1;

      } 

      if (shorter_linesegment_second_point_coords[1] > longer_line_first_point_coords[1]) {

	y_points_to_compare[0] = 1;
	y_points_to_compare[1] = 0;

      } 

      else {

	return null_vector;

      }

    }


    if (y_is_increasing_along_path_shorter_linesegment == false && y_is_increasing_along_path_longer_line == true) {

      if (shorter_linesegment_second_point_coords[1] < longer_line_first_point_coords[1]) {

	y_points_to_compare[0] = 1;
	y_points_to_compare[1] = 0;

      } 

      if (shorter_linesegment_first_point_coords[1] > longer_line_second_point_coords[1]) {

	y_points_to_compare[0] = 0;
	y_points_to_compare[1] = 1;

      } 

      else {

	return null_vector;

      }

    }

    if (y_is_increasing_along_path_shorter_linesegment == false && y_is_increasing_along_path_longer_line == false) {

      if (shorter_linesegment_second_point_coords[1] < longer_line_second_point_coords[1]) {

        y_points_to_compare[0] = 1;
        y_points_to_compare[1] = 1;

      }

      if (shorter_linesegment_first_point_coords[1] > longer_line_first_point_coords[1]) {

        y_points_to_compare[0] = 0;
        y_points_to_compare[1] = 0;

      }

      else {

        return null_vector;

      }

    }


    // z-coordinates

    if (z_is_increasing_along_path_shorter_linesegment == true && z_is_increasing_along_path_longer_line == true) {

      if (shorter_linesegment_first_point_coords[2] < longer_line_first_point_coords[2]) {

	z_points_to_compare[0] = 0;
	z_points_to_compare[1] = 0;

      } 

      if (shorter_linesegment_second_point_coords[2] > longer_line_second_point_coords[2]) {

	z_points_to_compare[0] = 1;
	z_points_to_compare[1] = 1;

      } 

      else {

	return null_vector;

      }

    }

    if (z_is_increasing_along_path_shorter_linesegment == true && z_is_increasing_along_path_longer_line == false) {

      if (shorter_linesegment_first_point_coords[2] < longer_line_second_point_coords[2]) {

        z_points_to_compare[0] = 0;
        z_points_to_compare[1] = 1;

      }

      if (shorter_linesegment_second_point_coords[2] > longer_line_first_point_coords[2]) {

        z_points_to_compare[0] = 1;
        z_points_to_compare[1] = 0;

      }

      else {

        return null_vector;

      }

    }

    if (z_is_increasing_along_path_shorter_linesegment == false && z_is_increasing_along_path_longer_line == true) {

      if (shorter_linesegment_second_point_coords[2] < longer_line_first_point_coords[2]) {

        z_points_to_compare[0] = 1;
        z_points_to_compare[1] = 0;

      }

      if (shorter_linesegment_first_point_coords[2] > longer_line_second_point_coords[2]) {

        z_points_to_compare[0] = 0;
        z_points_to_compare[1] = 1;

      }

      else {

        return null_vector;

      }

    }

    if (z_is_increasing_along_path_shorter_linesegment == false && z_is_increasing_along_path_longer_line == false) {

      if (shorter_linesegment_second_point_coords[2] < longer_line_second_point_coords[2]) {

        z_points_to_compare[0] = 1;
        z_points_to_compare[1] = 1;

      }

      if (shorter_linesegment_first_point_coords[2] > longer_line_first_point_coords[2]) {

        z_points_to_compare[0] = 0;
        z_points_to_compare[1] = 0;

      }

      else {

        return null_vector;

      }

    }

    // If the function has made it this far, then check to ensure that the components of the three lists are the same.
    for (size_t list_check = 0; list_check < 2; list_check++) {

      if (x_points_to_compare[list_check] != y_points_to_compare[list_check])
	return null_vector;

      if (x_points_to_compare[list_check] != z_points_to_compare[list_check])
	return null_vector;

      if (y_points_to_compare[list_check] != z_points_to_compare[list_check])
	return null_vector;

    }

    // The points that you want to return are:
    // (1) The point at the index ('0' or '1') from the 'shorter_linesegment'.  This the direction that you want to extend the track in.
    // (2) The point NOT at the index ('0' or '1') from the 'longer_line'.  This is the point on the 'longer_line' opposite the side that is being extended.

    // Declare the vector that will contain the output coordinates.
    std::vector < std::vector< double > > output_coordinate_vector;
    
    // Do this with a series of conditional statements.
    if (x_points_to_compare[0] == 0)
      output_coordinate_vector.push_back(shorter_linesegment_first_point_coords);

    if (x_points_to_compare[0] == 1) 
      output_coordinate_vector.push_back(shorter_linesegment_second_point_coords);

    if (x_points_to_compare[1] == 0) 
      output_coordinate_vector.push_back(longer_line_second_point_coords);

    if (x_points_to_compare[1] == 1)
      output_coordinate_vector.push_back(longer_line_first_point_coords);


    // Return 'output_coordinate_vector'
    return output_coordinate_vector;

  }

  // Define a function that will take the two coordinates of the points, find out which coordinates that they are on the short segment and long line endpoints, and adjust the correct track in the ‘track_v’ vector.

  // The inputs to the function are the following:   (1) ‘tracks_v’ - it should take in the vector of tracks that are being considered as being removed.
  //                                                 (2) ‘itrack_a’ - this is the index of ‘track_a’ in the ‘tracks_v’ vector.
  //                                                 (3) ’itrack_b’ - this is the index of ‘track_b’ in the ‘tracks_v’ vector.
  //                                                 (4) ’track_a’  - this is the ‘track_a’ object from the outer loop of the ‘process’ function.
  //                                                 (5) ‘track_b’   - this is the ‘track_b’ object form the inner loop of the ‘process’ function.
  //                                                 (6) ‘endpoint_coordinates_of_new_track’ - these are the endpoints of the new track that I will use to make the new track object.
  //  This function will be void, because it will just operate on the ‘tracks_v’ vector.
  void  concatenatingTwoTrackParts(std::vector< BMTrackCluster3D >& tracks_v, int itrack_a, int itrack_b, BMTrackCluster3D& track_a, BMTrackCluster3D& track_b, std::vector < std::vector <double> > endpoint_coordinates_of_new_track) {

    // Define vectors of doubles for the two sets of coordinates that are output by 'endpoint_coordinates_of_new_track'.
    // 'endpoint_coordinates_of_new_track[0]' corresponds to the 'shorter_lineseg' points on the extended track and ''endpoint_coordinates_of_new_track[1]' corresponds to the 'longer_line' points on the extended track.

    // 'push_back' the first entry of each of the vectors into this list, because those vectors contain vectors themselves (the coordinates of the two endpoints).
    std::vector <double> shorter_lineseg_endpoint_coords;
    shorter_lineseg_endpoint_coords.push_back(endpoint_coordinates_of_new_track[0][0]);

    std::vector <double> longer_line_endpoint_coords;
    longer_line_endpoint_coords.push_back(endpoint_coordinates_of_new_track[1][0]);

    // Find out which endpoint of the two tracks that this corresponds to.  The way this is done is by comparing the coordinates of one of the track endpoints with the coordinates of another track and its endpoints. 
    // I will use a system of indices to determine which track & which point these coordinates correspond to.  
    // First Index: 0 - 'track_a' coordinate, 1 - 'track_b' coordinate
    // Second Index: 0 - starting track coordinate, 1 - final track coordinate
    
    // I will do this by using a series of 'if' statements within a loop comparing the  x-coordinates of both the 'longer_line' and the 'shorter_linesegment' to the starting and ending x-coordinates of 'track_a' and 'track_b'.
    // Declare a vector of 'int's for each set of indices in this naming scheme.
    std::vector <int> shorter_lineseg_indices;
    std::vector <int> longer_line_indices;
          
    
    // I won't use equality to do this, but use a very small difference between the two coordinates.
    // Declare a vector that contains the two x-coordinates of 'longer_line' and 'shorter_lineseg'
    std::vector <double> endpoint_x_coordinates;
    endpoint_x_coordinates.push_back(shorter_lineseg_endpoint_coords[0]);
    endpoint_x_coordinates.push_back(longer_line_endpoint_coords[0]);
    
    // Declare 'int' variables for the indices that I will be using.
    int first_index  = 0;
    int second_index = 0;
          
    // Use a loop to assign the proper indices to the 'endpoint_x_coordinates'
    for (size_t i = 0; i < endpoint_x_coordinates.size(); i++) {

      // I am going to end up using nested 'if' loops, but it will not be terribly extensive.
      if (fabs(track_a.path3d.front()[0] - endpoint_x_coordinates[i]) < 0.001) {
	
	first_index  = 0;
	second_index = 0;

      }
      
      if (fabs(track_a.path3d.back()[0] - endpoint_x_coordinates[i]) < 0.001) {

	first_index  = 0;
	second_index = 1;
	
      }

      if (fabs(track_b.path3d.front()[0] - endpoint_x_coordinates[i]) < 0.001) {
	
	first_index  = 1;
	second_index = 0;

      }

      if (fabs(track_b.path3d.back()[0] - endpoint_x_coordinates[i]) < 0.001) {

	first_index  = 1;
	second_index = 1;

      }

      // Now, based on which index that we're using in the loop, fill either 'shorter_lineseg_indices' or 'longer_line_indices'
      if (i == 0) {

	shorter_lineseg_indices[0] = first_index;
	shorter_lineseg_indices[1] = second_index;

      }

      if (i == 1) {

	longer_line_indices[0] = first_index;
	longer_line_indices[1] = second_index;
	  
      }
      
    }

    
    // Declare information for all of the attributes of an object of type 'BMTrackCluster3D' that will have to appended onto the track belonging to the line.  Depending on the values for the indices that we find, these will be added on accordingly.

    // Info for 'short_lineseg'
    int shorter_lineseg_row;
    int shorter_lineseg_tick;
    std::vector <int> shorter_lineseg_wire;
    std::vector <double> shorter_lineseg_3Dpos;
    BoundarySpacePoint shorter_lineseg_endpts;
    BoundaryEnd_t shorter_lineseg_type;
    std::vector < BMTrackCluster2D > shorter_lineseg_plane_paths;
    std::vector< std::vector<double> > shorter_lineseg_path3d;
    int shorter_lineseg_track2d_index;

    // Info for 'longer_line'
    int longer_line_row;
    int longer_line_tick;
    std::vector<int> longer_line_wire;
    std::vector<double> longer_line_3Dpos;
    BoundarySpacePoint longer_line_endpts;
    BoundaryEnd_t longer_line_type;
    std::vector < BMTrackCluster2D > longer_line_plane_paths;
    std::vector< std::vector<double> > longer_line_path3d;
    int longer_line_track2d_index;
    
          

    // With this information, I can now tell which geometrical object, 'longer_line' or 'shorter_lineseg', corresponds to 'track_a' and which corresponds to 'track_b'.
    // The outer loop will consist of which track that the 'longer_line' belongs to.  That is the track that we will append information onto.
    // First, the case in which 'longer_line' corresponds to 'track_a'.
    if (longer_line_indices[0] == 0) {

      // 'longer_line' corresponds to 'track_a' and 'shorter_lineseg' corresponds to 'track_b'.
      
      // The three pieces of information that are specific to the track in general and not one of its endpoints, 'plane_paths', 'path3D', and 'track2D_index', can be set right now before we consider which point the endpoint corresponds to.

      // 'longer_line'
      longer_line_plane_paths         = track_a.plane_paths;
      longer_line_path3d              = track_a.path3d;
      longer_line_track2d_index       = track_a.track2d_index;
      
      // 'shorter_lineseg'
      shorter_lineseg_plane_paths     = track_b.plane_paths;
      shorter_lineseg_path3d          = track_b.path3d;
      shorter_lineseg_track2d_index   = track_b.track2d_index;
      
      
      // Within this loop, find out the coordinate of each of the tracks that the coordinates in each of the 'endpoint_coords' vectors correspond to.  There are two possibilities for each coordinate, and the orientation-specific (point-specific) information about each of the geometrical objects can be assigned in these loops.
    // First, that the endpoint on the 'longer_line' corresponds to the 'starting' point on 'track_a'.
      if (longer_line_indices[1] == 0) {
	     
	longer_line_row    = track_a.row_start;
	longer_line_tick   = track_a.tick_start;
	longer_line_wire   = track_a.start_wire;
	longer_line_3Dpos  = track_a.start3D;
	longer_line_endpts = track_a.start_endpts;
	longer_line_type   = track_a.start_type;

      }

      
      // Second, that the endpoint on the 'longer_line' corresponds to the 'ending' point on 'track_a'.
      if (longer_line_indices[1] == 1) {
	
	longer_line_row    = track_a.row_end;
	longer_line_tick   = track_a.tick_end;
	longer_line_wire   = track_a.end_wire;
	longer_line_3Dpos  = track_a.end3D;
	longer_line_endpts = track_a.end_endpts;
	longer_line_type   = track_a.end_type;
	
      }
      
      // Third, that the endpoint on the 'shorter_lineseg' corresponds to the 'starting' point on 'track_b'.
      if (shorter_lineseg_indices[1] == 0) {

	shorter_lineseg_row    = track_b.row_start;
	shorter_lineseg_tick   = track_b.tick_start;
	shorter_lineseg_wire   = track_b.start_wire;
	shorter_lineseg_3Dpos  = track_b.start3D;
	shorter_lineseg_endpts = track_b.start_endpts;
	shorter_lineseg_type   = track_b.start_type;
	
      }
      
      // Fourth and lastly, that the endpoint on the 'shorter_lineseg' corresponds to the 'ending' point on 'track_b'
      if (shorter_lineseg_indices[1] == 1) {

	shorter_lineseg_row    = track_b.row_end;
	shorter_lineseg_tick   = track_b.tick_end;
	shorter_lineseg_wire   = track_b.end_wire;
	shorter_lineseg_3Dpos  = track_b.end3D;
	shorter_lineseg_endpts = track_b.end_endpts;
	shorter_lineseg_type   = track_b.end_type;
	
      }
      
      
      // Now, based on the orientation of the y-coordinates of the two points, I can update the information in 'tracks_v.at(itrack_a)'.
      // I will make this track downwards going in the output.
      
      // First, the case in which the start of the track is the 'longer_line_endpoint_coordinates' and the end of the track corresponds to the 'shorter_lineseg_endpoint_coords'.
      if (longer_line_endpoint_coords[1] > shorter_lineseg_endpoint_coords[1]) {
	
	// Update the 'start' info first for the track with the 'longer_line' info.
	tracks_v.at(itrack_a).row_start    = longer_line_row;
	tracks_v.at(itrack_a).tick_start   = longer_line_tick;
	tracks_v.at(itrack_a).start_wire   = longer_line_wire;
	tracks_v.at(itrack_a).start3D      = longer_line_3Dpos;
	tracks_v.at(itrack_a).start_endpts = longer_line_endpts;
	tracks_v.at(itrack_a).start_type   = longer_line_type;
	
	
	// Update the 'end' info next for the track with the 'shorter_lineseg' info.
	tracks_v.at(itrack_a).row_end      = shorter_lineseg_row;
	tracks_v.at(itrack_a).tick_end     = shorter_lineseg_tick;
	tracks_v.at(itrack_a).end_wire     = shorter_lineseg_wire;
	tracks_v.at(itrack_a).end3D        = shorter_lineseg_3Dpos;
	tracks_v.at(itrack_a).end_endpts   = shorter_lineseg_endpts;
	tracks_v.at(itrack_a).end_type     = shorter_lineseg_type;

	// For the information for the entire track, you have to be a little careful.
	
	// The track index will just be the index of 'track_a'.  As a sanity check and for completeness, I will assign this again.
	tracks_v.at(itrack_a).track2d_index = longer_line_track2d_index;
	
	// Set 'tracks_v.at(itrack_a).plane_paths' equal to 'longer_line_plane_paths' just as a sanity check.
	tracks_v.at(itrack_a).plane_paths   = longer_line_plane_paths;
	
	// Set 'tracks_v.at(itrack_a).path3d' equal to 'longer_line_path3d' just as a sanity check.
	tracks_v.at(itrack_a).path3d        = longer_line_path3d;
	
	// Append the coordinates of 'shorter_lineseg_plane_paths' one at a time to 'tracks_v.at(itrack_a).plane_paths' AT THE END, because the 'shorter_lineseg' is at a lesser y-coordinate than 'longer_line' is.
	for (size_t plane_paths_iter = 0; plane_paths_iter < shorter_lineseg_plane_paths.size(); plane_paths_iter++) {
	    
	  tracks_v.at(itrack_a).plane_paths.push_back(shorter_lineseg_plane_paths[plane_paths_iter]);
	    
	}

	// Append the coordinates of 'shorter_lineseg_path3d' one at a time to 'tracks_v.at(itrack_a).path3D' AT THE END, because the 'shorter_lineseg' is a lesser y-coordinate than the 'longer_line' is.
	for (size_t path3d_iter = 0; path3d_iter < shorter_lineseg_path3d.size(); path3d_iter++) {
	    
	  tracks_v.at(itrack_a).path3d.push_back(shorter_lineseg_path3d[path3d_iter]);
	    
	}
	

      } // End the case in which 'longer_line' has a greater y-coordinate than 'shorter_lineseg'.

      // Next, the case in which the start of the track is the 'shorter_lineseg_endpoint_coordinates' and the end of the track corresponds to the 'longer_line_endpoint_coords'.
      if (shorter_lineseg_endpoint_coords[1] > longer_line_endpoint_coords[1]) {
	
	// Update the 'start' info first for the track with the 'shorter_lineseg' info.
	tracks_v.at(itrack_a).row_start    = shorter_lineseg_row;
	tracks_v.at(itrack_a).tick_start   = shorter_lineseg_tick;
	tracks_v.at(itrack_a).start_wire   = shorter_lineseg_wire;
	tracks_v.at(itrack_a).start3D      = shorter_lineseg_3Dpos;
	tracks_v.at(itrack_a).start_endpts = shorter_lineseg_endpts;
	tracks_v.at(itrack_a).start_type   = shorter_lineseg_type;

	// Update the 'end' info for the track with the 'longer_line' info.
	tracks_v.at(itrack_a).row_end      = longer_line_row;
	tracks_v.at(itrack_a).tick_end     = longer_line_tick;
	tracks_v.at(itrack_a).end_wire     = longer_line_wire;
	tracks_v.at(itrack_a).end3D        = longer_line_3Dpos;
	tracks_v.at(itrack_a).end_endpts   = longer_line_endpts;
	tracks_v.at(itrack_a).end_type     = longer_line_type;
	
	// For the information for the entire track, you have to be a little careful.
	
	// The track index will just be the index of 'track_a'.  As a sanity check and for completeness, I will assign this again.  
	tracks_v.at(itrack_a).track2d_index = longer_line_track2d_index;
	
	// Try this a different way....clear() the 'plane_paths' and 'path3d' vectors.
	tracks_v.at(itrack_a).plane_paths.clear();
	tracks_v.at(itrack_a).path3d.clear();

	// I will now fill these two vectors in the correct order.
	// Append the 'shorter_lineseg' information onto the front of both of the vectors and then append the 'longer_line' information onto the back of both of the vectors.
	
	// 'shorter_lineseg' 'plane_paths' information
	for (size_t shorter_lineseg_plane_paths_iter_first = 0; shorter_lineseg_plane_paths_iter_first < shorter_lineseg_plane_paths.size(); shorter_lineseg_plane_paths_iter_first++) {
	    
	  tracks_v.at(itrack_a).plane_paths.push_back(shorter_lineseg_plane_paths[shorter_lineseg_plane_paths_iter_first]);

	}

	// 'longer_line' 'plane_paths' information
	for (size_t longer_line_plane_paths_iter_first = 0; longer_line_plane_paths_iter_first < longer_line_plane_paths.size(); longer_line_plane_paths_iter_first++) {

	  tracks_v.at(itrack_a).plane_paths.push_back(longer_line_plane_paths[longer_line_plane_paths_iter_first]);

	}

	// 'shorter_lineseg' 'path3d' information
	for (size_t shorter_lineseg_path3d_iter_first = 0; shorter_lineseg_path3d_iter_first < shorter_lineseg_path3d.size(); shorter_lineseg_path3d_iter_first++) {
	    
	  tracks_v.at(itrack_a).path3d.push_back(shorter_lineseg_path3d[shorter_lineseg_path3d_iter_first]);

	}

	// 'longer_line' 'path3d' information
	for (size_t longer_line_path3d_iter_first = 0; longer_line_path3d_iter_first < longer_line_path3d.size(); longer_line_path3d_iter_first++) {

	  tracks_v.at(itrack_a).path3d.push_back(longer_line_path3d[longer_line_path3d_iter_first]);

	} 

      } // End the case in which the 'shorter_lineseg' is at a greater y-coordinate than the 'longer_line'
            
    } // End the case in which 'longer_line' corresponds to 'track_a'.

    
    // Second, the case in which 'longer_line' corresponds to 'track_b'.                                
    if (longer_line_indices[0] == 1) {

      // 'shorter_lineseg' corresponds to 'track_a' and 'longer_line' corresponds to 'track_b'.                                               
      
      // The three pieces of information that are specific to the track in general and not one of its endpoints, 'plane_paths', 'path3D', 
      // and 'track2d_index', can be set right now before we consider which point the endpoint corresponds to.                           

      // 'shorter_lineseg'
      shorter_lineseg_plane_paths         = track_a.plane_paths;
      shorter_lineseg_path3d              = track_a.path3d;
      shorter_lineseg_track2d_index       = track_a.track2d_index;

      // 'longer_line'
      longer_line_plane_paths             = track_b.plane_paths;
      longer_line_path3d                  = track_b.path3d;
      longer_line_track2d_index           = track_b.track2d_index;

      // Within this loop, find out the coordinate of each of the tracks that the coordinates in each of the 'endpoint_coords' vectors correspond to.  There are two possibilities for each coordinate, and the orientation-specific (point-specific) information about each of the geometrical objects can be assigned in these loops.                                                                                                                                                                                  
      // First, that the endpoint on the 'shorter_lineseg' corresponds to the 'starting' point on 'track_a'. 
      if (shorter_lineseg_indices[1] == 0) {

	shorter_lineseg_row    = track_a.row_start;
    	shorter_lineseg_tick   = track_a.tick_start;
     	shorter_lineseg_wire   = track_a.start_wire;
    	shorter_lineseg_3Dpos  = track_a.start3D;
     	shorter_lineseg_endpts = track_a.start_endpts;
     	shorter_lineseg_type   = track_a.start_type;
	
      }
      
      // Second, that the endpoint on the 'shorter_lineseg' corresponds to the 'ending' point on 'track_a'
      if (shorter_lineseg_indices[1] == 1) {

	shorter_lineseg_row    = track_a.row_end;
    	shorter_lineseg_tick   = track_a.tick_end;
     	shorter_lineseg_wire   = track_a.end_wire;
     	shorter_lineseg_3Dpos  = track_a.end3D;
     	shorter_lineseg_endpts = track_a.end_endpts;
     	shorter_lineseg_type   = track_a.end_type;

      }

      // Third, that the endpoint on the 'longer_line' corresponds to the 'starting' point on 'track_b'
      if (longer_line_indices[1] == 0) {
	  
    	longer_line_row        = track_b.row_start;
    	longer_line_tick       = track_b.tick_start;
     	longer_line_wire       = track_b.start_wire;
    	longer_line_3Dpos      = track_b.start3D;
     	longer_line_endpts     = track_b.start_endpts;
     	longer_line_type       = track_b.start_type;
	
      }

      // Fourth and lastly, that the endpoint on the 'longer_line' corresponds to the 'ending' point on 'track_b'
      if (longer_line_indices[1] == 1) {
	
     	longer_line_row        = track_b.row_end;
     	longer_line_tick       = track_b.tick_end;
     	longer_line_wire       = track_b.end_wire;
     	longer_line_3Dpos      = track_b.end3D;
     	longer_line_endpts     = track_b.end_endpts;
     	longer_line_type       = track_b.end_type;
	
      }

      // Now, based on the orientation of the y-coordinates of the two points, I can update the information in 'tracks_v.at(itrack_b)'. 
      // I will make this track downwards going in the output.                                                                             
      
      // First, the case in which the start of the track is the 'longer_line_endpoint_coordinates' and the end of the track corresponds 
      // to the 'shorter_lineseg_endpoint_coords'.                                                                                                                                                      
      if (longer_line_endpoint_coords[1] > shorter_lineseg_endpoint_coords[1]) {
	  
     	// Update the 'start' info first for the track with the 'longer_line' info.
     	tracks_v.at(itrack_b).row_start    = longer_line_row;
     	tracks_v.at(itrack_b).tick_start   = longer_line_tick;
    	tracks_v.at(itrack_b).start_wire   = longer_line_wire;
     	tracks_v.at(itrack_b).start3D      = longer_line_3Dpos;
    	tracks_v.at(itrack_b).start_endpts = longer_line_endpts;
     	tracks_v.at(itrack_b).start_type   = longer_line_type;
	
     	// Update the 'end' info second for the track with the 'shorter_lineseg' info.
     	tracks_v.at(itrack_b).row_end      = shorter_lineseg_row;
     	tracks_v.at(itrack_b).tick_end     = shorter_lineseg_tick;
     	tracks_v.at(itrack_b).end_wire     = shorter_lineseg_wire;
    	tracks_v.at(itrack_b).end3D        = shorter_lineseg_3Dpos;
     	tracks_v.at(itrack_b).end_endpts   = shorter_lineseg_endpts;
    	tracks_v.at(itrack_b).end_type     = shorter_lineseg_type;
	
     	// Set the information general to the entire track now.
	
     	// As a sanity check, set 'tracks_v.at(itrack_b).track2d_index' equal to 'longer_line_index' (it should already have that value, but this is for completeness).
    	tracks_v.at(itrack_b).track2d_index = longer_line_track2d_index;

     	// As another sanity check, set 'tracks_v.at(itrack_b).plane_paths' equal to 'longer_line_plane_paths'
     	tracks_v.at(itrack_b).plane_paths   = longer_line_plane_paths;
	
    	// As a final sanity check, set 'tracks_v.at(itrack_b).path3d' equal to 'longer_line_path3d'
     	tracks_v.at(itrack_b).path3d        = longer_line_path3d;
	
     	// Append information from 'shorter_lineseg_plane_paths' onto the end of the 'tracks_v.at(itrack_b).plane_paths' vector, starting with the first entry.
     	for (size_t plane_paths_iter = 0; plane_paths_iter < shorter_lineseg_plane_paths.size(); plane_paths_iter++) {
          
     	  tracks_v.at(itrack_b).plane_paths.push_back(shorter_lineseg_plane_paths[plane_paths_iter]);
	    
    	}
	
     	// Append information from 'shorter_line_seg_path3d' onto the end of the 'tracks_v.at(itrack_b).path3d' vector, starting with the first entry.
	
     	for (size_t path3d_iter = 0; path3d_iter < shorter_lineseg_path3d.size(); path3d_iter++) {
	    
    	  tracks_v.at(itrack_b).path3d.push_back(shorter_lineseg_path3d[path3d_iter]);
	    
	}
	
      }

      // Second, the case in which the start of the track is at 'shorter_lineseg_endpoint_coords' and the end of the track is located at 'longer_line_endpoint_coords'.
      if (shorter_lineseg_endpoint_coords[1] > longer_line_endpoint_coords[1]) {
	
   	// Update the 'start' info first for the track with the 'shorter_lineseg' info.
   	tracks_v.at(itrack_b).row_start    = shorter_lineseg_row;
     	tracks_v.at(itrack_b).tick_start   = shorter_lineseg_tick;
     	tracks_v.at(itrack_b).start_wire   = shorter_lineseg_wire;
     	tracks_v.at(itrack_b).start3D      = shorter_lineseg_3Dpos;
     	tracks_v.at(itrack_b).start_endpts = shorter_lineseg_endpts;
     	tracks_v.at(itrack_b).start_type   = shorter_lineseg_type;
	
	// Update the 'end' info second for the track with the 'longer_line' info.
	tracks_v.at(itrack_b).row_end      = longer_line_row;
	tracks_v.at(itrack_b).tick_end     = longer_line_tick;
	tracks_v.at(itrack_b).end_wire     = longer_line_wire;
	tracks_v.at(itrack_b).end3D        = longer_line_3Dpos;
	tracks_v.at(itrack_b).end_endpts   = longer_line_endpts;
	tracks_v.at(itrack_b).end_type     = longer_line_type;
	
	// Set the information general to the entire track now.
	
	// As a sanity check, set 'tracks_v.at(itrack_b).track2d_index' equal to 'longer_line_index' (it should 
	// already have that value, but this is for completeness).                                                                          
	tracks_v.at(itrack_b).track2d_index = longer_line_track2d_index;
	
	// Clear out the 'plane_paths' and 'path3D' components of 'tracks_v.at(itrack_b)'.
	tracks_v.at(itrack_b).plane_paths.clear();
	tracks_v.at(itrack_b).path3d.clear();
	
	// Append information from 'shorter_lineseg_plane_paths' onto the 'plane_paths' vector FIRST.
	for (size_t shorter_lineseg_plane_paths_iter_second = 0; shorter_lineseg_plane_paths_iter_second < shorter_lineseg_plane_paths.size(); shorter_lineseg_plane_paths_iter_second++) {
	    
	  tracks_v.at(itrack_b).plane_paths.push_back(shorter_lineseg_plane_paths[shorter_lineseg_plane_paths_iter_second]);
	    
	}

	// Append information from 'longer_line_plane_paths' onto the 'plane_paths' vector SECOND
	for (size_t longer_line_plane_paths_iter_second = 0; longer_line_plane_paths_iter_second < longer_line_plane_paths.size(); longer_line_plane_paths_iter_second++) {

	  tracks_v.at(itrack_b).plane_paths.push_back(longer_line_plane_paths[longer_line_plane_paths_iter_second]);

	}

	// Append information from 'shorter_lineseg_path3d' onto the 'path3d' vector FIRST.
	for (size_t shorter_lineseg_path3d_iter_second = 0; shorter_lineseg_path3d_iter_second < shorter_lineseg_path3d.size(); shorter_lineseg_path3d_iter_second++) {

	  tracks_v.at(itrack_b).path3d.push_back(shorter_lineseg_path3d[shorter_lineseg_path3d_iter_second]);

	}

	// Append information from 'longer_line_path3d' onto the 'path3d' vector SECOND.
	for (size_t longer_line_path3d_iter_second = 0; longer_line_path3d_iter_second < longer_line_path3d.size(); longer_line_path3d_iter_second++) {

	  tracks_v.at(itrack_b).path3d.push_back(longer_line_path3d[longer_line_path3d_iter_second]);

	}
	    
      } // End of the case in which the 'shorter_lineseg' has a greater y-coordinate than the 'longer_line'.
	
    } // End of the case in which the 'longer_line' (the point in the vector onto which items are appended) corresponds to 'track_b' in the loop.
      
  } // End of the function 'concatenatingTwoTrackParts'.


  // Generate the points on both the linesegment and the line.  The number of points will be equal to resolution value.                                                         
  // Inputs:  (1) 'point_one_coords' - The coordinates of one endpoint on the track.                                                                                            
  //          (2) 'point_two_coords' - The coordinates of another endpoint on the track.                                                                                        
  //          (3) 'shorter_dim_idx_vector' - The index in the coordinate vector of the two dimensions that are not greatest in spread.                                         
  //          (4) 'longest_dim_idx' - The index in the coordinate vector of the track dimension with the greatest spread.                                                      
  std::vector < std::vector<double> > Linear3DPostProcessor::getTrackPoints(const std::vector<double> point_one_coords, const std::vector<double> point_two_coords, int resolution, std::vector <int> shorter_dim_idx_vector, int longest_dim_idx) {

    // Divide the coordinates out into the ones corresponding to the dimension with the largest spread in each of the lines (the independent variable) and the ones corresponding to less of a spread in each of the lines (the dependent variables).                                                                                                              
    // Right now, I do not know which dimension in 3D Cartesian space these points correspond to.                                                                                   

    // First set of points                                                                                                                                                      
    double largest_spread_var_coord_first  = point_one_coords[longest_dim_idx];
    double shorter_spread_coord1_var_first = point_one_coords[shorter_dim_idx_vector[0]];
    double shorter_spread_coord2_var_first = point_one_coords[shorter_dim_idx_vector[1]];

    // Second set of points                                                                                                                                                 
    double largest_spread_var_coord_second   = point_two_coords[longest_dim_idx];
    double shorter_spread_coord1_var_second  = point_two_coords[shorter_dim_idx_vector[0]];
    double shorter_spread_coord2_var_second  = point_two_coords[shorter_dim_idx_vector[1]];

    // Divide the difference between 'largest_spread_var_first' and 'largest_spread_var_second' into 'resolution' number of equal segments, and then find which of the 'largest_spre\ad' variables is greatest.                                                                                                                                                 
    double segment_length = fabs(largest_spread_var_coord_first - largest_spread_var_coord_second);

    double lesser_largest_spread_coord  = largest_spread_var_coord_first;
    double greater_largest_spread_coord = largest_spread_var_coord_second;

    // If 'largest_spread_var_coord_first' is greater than 'largest_spread_var_coord_second', then you can redefine these variables according to their different relative magnitudes.                                                                                                                                                                               
    // This is important in moving from one point on the line to the other, but NOT in calculating the slope (the calculation of the slope is blind to the relative sign in the coordinates in the numerator and the denominator, because that cancels out).  This does not correlate between the coordinates on the line as well.                                     
    if (largest_spread_var_coord_first > largest_spread_var_coord_second) {

      lesser_largest_spread_coord  = largest_spread_var_coord_second;
      greater_largest_spread_coord = largest_spread_var_coord_first;

    }

    // Define the slopes and the y-intercepts for the two lines that you will be used.                                                                                         
    // Note: This is the only point at which I will use the coordinates in the 'shorter_spread_coord' variables, because those will be the outputs of these two lines.          
    double slope_line1 = (shorter_spread_coord1_var_first - shorter_spread_coord1_var_second)/(largest_spread_var_coord_first - largest_spread_var_coord_second);

    double slope_line2 = (shorter_spread_coord2_var_first - shorter_spread_coord2_var_second)/(largest_spread_var_coord_first- largest_spread_var_coord_second);

    double yint_line1 = shorter_spread_coord1_var_first - (slope_line1 * largest_spread_var_coord_first);

    double yint_line2 = shorter_spread_coord2_var_first - (slope_line2 * largest_spread_var_coord_first);

    // Next, I will divide the difference between 'greater_largest_spread_pt' and 'lesser_largest_spread_pt' into 'resolution' number of segments.                           
    double largest_spread_segment_length = (greater_largest_spread_coord - lesser_largest_spread_coord) / resolution;

    // Declare a vector of type 'double' for each coordinate that will contain 'resolution + 1' entries each.                                                            
    // I could also declare these vectors as empty and then append each of the values of these coordinates, but I will give them a finite size to make it easier to debug     
    // the algorithm if I run into problems.                                                                                                                                   
    std::vector <double> largest_spread_coord_values  (resolution + 1, 0.0);
    std::vector <double> shorter_spread_coord1_values (resolution + 1, 0.0);
    std::vector <double> shorter_spread_coord2_values (resolution + 1, 0.0);

    // Set a variable for the current value of 'largest_spread_coord', 'current_val_largest_spread_coord', and initialize it to 'lesser_largest_spread_coord'.                
    double current_val_largest_spread_coord = lesser_largest_spread_coord;

    // Initialize variables for the current value of 'shorter_spread_coord1' and 'shorter_spread_coord2', and initialize them both to 0.                                       
    double current_val_shorter_spread_coord1 = 0.;
    double current_val_shorter_spread_coord2 = 0.;

    // A shortcoming of this method: I may not exactly hit the last point on the line that the function is looping over right now.  That is because, with doubles, continuously adding increments of 'largest_spread_segment_length' may not exactly equal 'greater_largest_spread_coord' in the end.  In principle it should, but programming with doubles means that I may miss the final value by a fraction of a cm (which will not change the result in the end).                                                                                 
    // Begin a 'for' loop to fill the vectors declared above.                                                                                                               
    for (size_t line_iter = 0; line_iter < resolution + 1; line_iter++) {

      // Calculate 'current_val_shorter_spread_coord1' and 'current_val_shorter_spread_coord2'                                                                           
      // 'current_val_shorter_spread_coord1'                                                                                           
      current_val_shorter_spread_coord1 = (slope_line1*current_val_largest_spread_coord) + yint_line1;

      // 'current_val_shorter_spread_coord2'                                                                                                                               
      current_val_shorter_spread_coord2 = (slope_line2*current_val_largest_spread_coord) + yint_line2;
      
      // Set the value of the respective vector for each of the three variables at index number 'line_iter' to the current_value of the variable.                          
      largest_spread_coord_values[line_iter]  = current_val_largest_spread_coord;
      shorter_spread_coord1_values[line_iter] = current_val_shorter_spread_coord1;
      shorter_spread_coord2_values[line_iter] = current_val_shorter_spread_coord2;

    }

    // After the three vectors have been filled, I will use a loop to fill the final vector, 'line_coord_values', which will be an outer vector of inner vectors, with the inner vectors being the [x, y, z] coordinates at each particular point.  To do this, I must know which index 'longest_dim_idx' and 'shorter_dim_idx_vector' correspond to.            

    // Fill this vector slightly differently for the time being.                                                                                                               
    std::vector < std::vector <double> > line_coord_values;

    // Next, begin a loop of length 'largest_spread_coord_values' (which will serve as a test that the other vectors are the same length as well.                           
    for (size_t coord_iter = 0; coord_iter < largest_spread_coord_values.size(); coord_iter++) {

      // 'push_back' another element onto 'line_coord_values'
      line_coord_values.push_back(std::vector< double > (3, 0.0));

      // Begin the cases for which coordinates 'longest_dim_idx' and the 'shorter_dim_idx_vector' correspond to.                                                          
      // longest dimension - x, shorter dimension #1 - y, shorter dimension #2 - z                                                                                         
      if (longest_dim_idx == 0 && shorter_dim_idx_vector[0] == 1 && shorter_dim_idx_vector[1] == 2) {

	// Set the information in the last entry of the new vector.
	line_coord_values[line_coord_values.size() - 1][0] = largest_spread_coord_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = shorter_spread_coord1_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = shorter_spread_coord2_values[coord_iter]; 
      
      }

      // longest dimension - x, shorter dimension #1 - z, shorter dimension #2 - y                                                                                          
      else if (longest_dim_idx == 0 && shorter_dim_idx_vector[0] == 2 && shorter_dim_idx_vector[1] == 1) {

	// Set the information in the last entry of the new vector.
	line_coord_values[line_coord_values.size() - 1][0] = largest_spread_coord_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = shorter_spread_coord2_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = shorter_spread_coord1_values[coord_iter];

      }

      // longest dimension - y, shorter dimension #1 - x, shorter dimension #2 - z                                                                                        
      else if (longest_dim_idx == 1 && shorter_dim_idx_vector[0] == 0 && shorter_dim_idx_vector[1] == 2) {

	line_coord_values[line_coord_values.size() - 1][0] = shorter_spread_coord1_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = largest_spread_coord_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = shorter_spread_coord2_values[coord_iter];

      }

      // longest dimension - y, shorter dimension #1 - z, shorter dimension #2 - x                                                                                          
      else if (longest_dim_idx == 1 && shorter_dim_idx_vector[0] == 2 && shorter_dim_idx_vector[1] == 0) {
	
	line_coord_values[line_coord_values.size() - 1][0] = shorter_spread_coord2_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = largest_spread_coord_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = shorter_spread_coord1_values[coord_iter];

      }

      // longest dimension - z, shorter dimension #1 - x, shorter dimension #2 - y                                                                                          
      else if (longest_dim_idx == 2 && shorter_dim_idx_vector[0] == 0 && shorter_dim_idx_vector[1] == 1) {

	line_coord_values[line_coord_values.size() - 1][0] = shorter_spread_coord1_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = shorter_spread_coord2_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = largest_spread_coord_values[coord_iter];

      }

      // longest dimension - z, shorter dimension #1 - y, shorter dimension #2 - x                                                                               
      else if (longest_dim_idx == 2 && shorter_dim_idx_vector[0] == 1 && shorter_dim_idx_vector[1] == 0) {

	line_coord_values[line_coord_values.size() - 1][0] = shorter_spread_coord2_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][1] = shorter_spread_coord1_values[coord_iter];
	line_coord_values[line_coord_values.size() - 1][2] = largest_spread_coord_values[coord_iter];

      }

    }

    // The points on the line have now been set.  I can return 'line_coord_values' and conclude the function.                                           
    return line_coord_values;

  }





  // Find out which one of the dimensions has the largest spread of the three dimensions.                                                                           
  // Input: (1) The Coordinates of the First Point on the track-like entity.                                                                                                
  //        (2) The Coordinates of the Second Point on the track-like entity.                                                                      
  int Linear3DPostProcessor::longestDimFinder(const std::vector<double> point_one_coords, const std::vector<double> point_two_coords) {

    // Define a variable for the smallest dimension.  Initialize it to the spread in the x-direction.                                                  
    double longest_dim     = fabs(point_one_coords[0] - point_two_coords[0]);

    // The index notation that I will use for the dimensions with respect to one another is '0' - x, '1' - y, '2' - z.                             
    int    longest_dim_idx = 0;

    // Use a 'for' loop to find which is the longest dimension.                                                                                                 
    for (int dim_i = 1; dim_i < 3; dim_i++) {

      // Check to see if the entry at this index has greater spread than the dimension at current value of 'longest_dim_idx'.                                      
      if (fabs(point_one_coords[dim_i] - point_two_coords[dim_i]) > longest_dim) {

        // Reset the value of 'longest_dim'.                                                                                                                            
	longest_dim     = fabs(point_one_coords[dim_i] - point_two_coords[dim_i]);

        // Reset the value of 'longest_dim_idx'.                                                                                                       
	longest_dim_idx = dim_i;

      }

    }

    // Once this loop is complete, you can return the integer value of 'longest_dim_idx'.                                                              
    return longest_dim_idx;

      }

}
