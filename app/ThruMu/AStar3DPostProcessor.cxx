/* Write an algorithm that will look through the points on the AStar3D nodes to see if they contain charge.
 *These paths will not necessarily be straight, which is why this algorithm differs from the 'Linear3DPostProcessor' module.

 * Author: Christopher Barnes (cbarnes4@fnal.gov)
 */

// Add an 'include' statement for the '.h' file in the header
#include "AStar3DPostProcessor.h"
#include "GeoAlgo.h"
#include "GeoLineSegment.h"

// Include the 'BMTrackCluster3D' source file so that I can access 'markImageWithTrack' instead of redefining it here.
#include "BMTrackCluster3D.h"


// Include 'math.h' so I will certainly have access to the 'sqrt' functionality that it has.
#include <math.h>

namespace larlitecv {

  // Define a function that will look at tracks that were generated with the 'AStar' algorithm, meaning that they are more curved, and try to separate them.
  // Input for this function: (1) 'tracks_v' - This the list of remaining tracks, of type 'BMTrackCluster', which contain the 
  std::vector< BMTrackCluster3D > AStar3DPostProcessor::separateAStarTracks(std::vector < BMTrackCluster3D >& tracks_v, double maximum_distance, const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs, const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size) {

    // Start a loop over the number of tracks in the 'tracks_v' array.
    const size_t n_tracks = tracks_v.size();

    // Introduce the vector of tracks that you want to skip because they have a maximum overlap with another track
    std::vector<bool> track_excluded( n_tracks, false );

    // Initialize the loop that same way that it is done in the 'Linear3DPostProcessor' algorithm.
    for (size_t itrack_a=0; itrack_a<n_tracks; itrack_a++) {

      // If this track has already been excluded in the inner loop, then you can continue.
      if ( track_excluded.at(itrack_a) ) continue;

      // Declare the 'BMTrackCluster3D' object at this point in the array.  
      BMTrackCluster3D& track_a = tracks_v.at(itrack_a);

      // We are concerned with this entire path, which we will trace out over a set of the pixels to see if they contain as much charge as this one does.
      for (size_t itrack_b=0; itrack_b<n_tracks; itrack_b++) {

	// If this track has already been excluded in the outer loop, then you can continue.
	if ( track_excluded.at(itrack_b) ) continue;

	// Immediately, if these two tracks are the same (meaning that they have the same index), then you can continue.
	if (itrack_a == itrack_b)
	  continue;

	// Declare the 'BMTrackCluster3D' object at this point in the array.
	BMTrackCluster3D& track_b = tracks_v.at(itrack_b);

	// Go through the track endpoints to see which endpoints are most extreme in each of the directions.

	// While these tracks by no means will be straight, I can assume that the relative placement of their starting points and their ending points to one another will indicate the direction of the track's progression within the TPC. 


	// Declare, for the track at the outer loop corresponding to index 'itrack_a', the 'starting_coordinates' as the coordinates least in their variable and the ending coordinates as the coordinates greatest in their variable.
	double least_x_value_on_track_end_a    = track_a.path3d.front()[0];
	double least_y_value_on_track_end_a    = track_a.path3d.front()[1];
	double least_z_value_on_track_end_a    = track_a.path3d.front()[2];

	double greatest_x_value_on_track_end_a = track_a.path3d.back()[0];
	double greatest_y_value_on_track_end_a = track_a.path3d.back()[1];
	double greatest_z_value_on_track_end_a = track_a.path3d.back()[2];

	// If these coordinates are not defined in the correct order, then you can redefine them in the right order.
	if (track_a.path3d.back()[0] < track_a.path3d.front()[0]) {

	  least_x_value_on_track_end_a    = track_a.path3d.back()[0];
	  greatest_x_value_on_track_end_a = track_a.path3d.front()[0];

	}

        if (track_a.path3d.back()[1] < track_a.path3d.front()[1]) {

          least_y_value_on_track_end_a    = track_a.path3d.back()[1];
          greatest_y_value_on_track_end_a = track_a.path3d.front()[1];

        }

        if (track_a.path3d.back()[2] < track_a.path3d.front()[2]) {

          least_z_value_on_track_end_a    = track_a.path3d.back()[2];
          greatest_z_value_on_track_end_a = track_a.path3d.front()[2];

        }

	// Repeat the same process for the 'track_b' coordinates
	double least_x_value_on_track_end_b    = track_b.path3d.front()[0];
        double least_y_value_on_track_end_b    = track_b.path3d.front()[1];
        double least_z_value_on_track_end_b    = track_b.path3d.front()[2];

        double greatest_x_value_on_track_end_b = track_b.path3d.back()[0];
        double greatest_y_value_on_track_end_b = track_b.path3d.back()[1];
        double greatest_z_value_on_track_end_b = track_b.path3d.back()[2];

        // If these coordinates are not defined in the correct order, then you can redefine them in the right order.                                                                 
        if (track_b.path3d.back()[0] < track_b.path3d.front()[0]) {

          least_x_value_on_track_end_b    = track_b.path3d.back()[0];
          greatest_x_value_on_track_end_b = track_b.path3d.front()[0];

        }

        if (track_b.path3d.back()[1] < track_b.path3d.front()[1]) {

          least_y_value_on_track_end_b    = track_b.path3d.back()[1];
          greatest_y_value_on_track_end_b = track_b.path3d.front()[1];

        }

        if (track_b.path3d.back()[2] < track_b.path3d.front()[2]) {

          least_z_value_on_track_end_b    = track_b.path3d.back()[2];
          greatest_z_value_on_track_end_b = track_b.path3d.front()[2];

        }

	// With correctly oriented coordinates, calculate the distance between one point and the other.
	double distance_between_lesser_value_points = sqrt((least_x_value_on_track_end_a - least_x_value_on_track_end_b)*(least_x_value_on_track_end_a - least_x_value_on_track_end_b) + (least_y_value_on_track_end_a - least_y_value_on_track_end_b)*(least_y_value_on_track_end_a - least_y_value_on_track_end_b) + (least_z_value_on_track_end_a - least_z_value_on_track_end_b)*(least_z_value_on_track_end_a - least_z_value_on_track_end_b));

	double distance_between_greater_value_points = sqrt((greatest_x_value_on_track_end_a - greatest_x_value_on_track_end_b)*(greatest_x_value_on_track_end_a - greatest_x_value_on_track_end_b) + (greatest_y_value_on_track_end_a - greatest_y_value_on_track_end_b)*(greatest_y_value_on_track_end_a - greatest_y_value_on_track_end_b) + (greatest_z_value_on_track_end_a - greatest_z_value_on_track_end_b)*(greatest_z_value_on_track_end_a - greatest_z_value_on_track_end_b));

	// Check to ensure that these two values are both less than the 'maximum_distance' defined above.
	// Continue if either of them are greater than this value.  That means that these are clearly distinct tracks.
	if ((distance_between_lesser_value_points > maximum_distance) || (distance_between_greater_value_points > maximum_distance)) {

	  continue;

	}

	// If you have made it this far, then mark up blank images that have charge above threshold for each of the tracks.
	

	// Here, I will set 'markedimgs' to [] so that it will be set to an empty image.  I will also 'markedvalue' in the 'markImageWithTrack' function as '1.0'.
	std::vector<larcv::Image2D> markedimgs_tracka;
	std::vector<larcv::Image2D> markedimgs_trackb;
	bool markvalue = 1.0;

	// Declare an object of class 'BMTrackCluster3D' in order to mark the two images created above.
	BMTrackCluster3D marker_obj;

	// 'track_a'
	// The 'marker_obj' can call the function 'markImageWithTrack', which is void but fills 'markedimgs_tracka'.
	marker_obj.markImageWithTrack(imgs, badchimgs, thresholds, neighborhood_size, markedimgs_tracka, markvalue);
	
	// 'track_b' (how can it tell which track that I'm referring to)
	// The 'marker_obj' can call the function 'markImageWithTrack', which is void but fills 'markedimgs_trackb'.
	marker_obj.markImageWithTrack(imgs, badchimgs, thresholds, neighborhood_size, markedimgs_trackb, markvalue);

	// Now, declare 'int' values for the number of points in each of the images that are above threshold

	// Declare variables for the number of pixels across each of the three planes that contain charge.
	int num_of_pixels_three_planes_tracka = 0;
	int num_of_pixels_three_planes_trackb = 0;

	
        // Both of these images are the same size, so I can count the number of points in the same loop.                                                                            
	// Loop over each of the plane images in the 'markedimgs_tracka' vector.
	for (size_t plane_iter = 0; plane_iter < markedimgs_tracka.size(); plane_iter++) {

	  // Loop over the rows and columns in the image.  Each of the plane images for both 'track_a' and 'track_b' should have the same number of rows and 
	  // columns, so I should not have vector capacity issues.
	  for (size_t row_iter = 0; row_iter < markedimgs_tracka[plane_iter].meta().rows(); row_iter++) {

	    for (size_t col_iter = 0; col_iter < markedimgs_tracka[plane_iter].meta().cols(); col_iter++) {

	      // Increment 'num_of_pixels_three_planes_tracka' if this entry is above threshold.
	      // I'll use a 'greater than' operator for '0.9' just to ensure that I do not run into errors with floating point equality comparisons.
	      if (markedimgs_tracka[plane_iter].pixel(row_iter, col_iter) > 0.9) {

		num_of_pixels_three_planes_tracka += 1;

	      }

	      // Increment 'num_of_pixels_three_planes_trackb' if this entry is above threshold.
	      if (markedimgs_trackb[plane_iter].pixel(row_iter, col_iter) > 0.9) {

		num_of_pixels_three_planes_trackb += 1;

	      }

	    }

	  }

	}

	// Compare the values of 'num_of_pixels_three_planes_tracka' and 'num_of_pixels_three_planes_trackb' to one another to see which image has more pixels that were equal to 1.0.
	if (num_of_pixels_three_planes_tracka > num_of_pixels_three_planes_trackb) {

	  // Set the value of 'tracks_excluded' at 'track_b' equal to 'true'
	  track_excluded[itrack_b] = true;

	  // Continue on to the next value of 'track_b' in the loop.
	  continue;

	}

	if (num_of_pixels_three_planes_trackb > num_of_pixels_three_planes_tracka) {
	  
	  // Set the value of 'tracks_excluded' at 'track_a' equal to 'true'
	  track_excluded[itrack_a] = true;

	  // Continue on to the next value of 'track_a' in the outer loop, meaning that you must break the inner loop.
	  break;

	}

      }

    }

    // collect the output                                                                                                                                                
    std::vector< BMTrackCluster3D > output_tracks;

    for ( size_t itrack=0; itrack<track_excluded.size(); itrack++ ) {
      if ( track_excluded.at(itrack) )
	continue;

      output_tracks.emplace_back( std::move(tracks_v.at(itrack)) );
    }

    tracks_v.clear();

    return output_tracks;

  }

}





			   
    


	  
