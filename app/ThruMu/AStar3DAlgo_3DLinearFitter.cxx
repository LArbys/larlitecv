#include "AStar3DAlgo.h"
#include <vector>
#include <cmath>

// larcv
#include "UBWireTool/UBWireTool.h"

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace larlitecv {

  // I have added the variables 'time_comp_factor' and 'wire_comp_factor' to this definition when they were not here previously
  std::vector<AStar3DNode> AStar3DAlgo_3DFitter::findpath( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const int start_row, const int goal_row, const std::vector<int>& start_cols, const std::vector<int>& goal_cols  ) {
    
    const larcv::ImageMeta& meta = img_v.front().meta();

    if ( verbose>0 ) 
      std::cout << "[[ASTAR 3D ALGO START]]" << std::endl;

    // turn image pixel coordinates into a 3D start and goal point
    std::vector<int> start_wids(start_cols.size());
    std::vector<int> goal_wids(goal_cols.size());    
    for (size_t p=0; p<goal_cols.size(); p++) {
      start_wids[p] = img_v.at(p).meta().pos_x( start_cols[p] );
      goal_wids[p]  = img_v.at(p).meta().pos_x( goal_cols[p] );      
    }
    double start_tri = 0.;
    int start_crosses = 0;
    std::vector< float > poszy_start(2,0.0);
    larcv::UBWireTool::wireIntersection( start_wids, poszy_start, start_tri, start_crosses );

    double goal_tri = 0.;
    int goal_crosses = 0;
    std::vector< float > poszy_goal(2,0.0);
    larcv::UBWireTool::wireIntersection( goal_wids, poszy_goal, goal_tri, goal_crosses );

    std::vector<float> startpos(3,0);
    startpos[0] = meta.pos_y( start_row );
    startpos[1] = poszy_start[1];
    startpos[2] = poszy_start[0];

    std::vector<float> goalpos(3,0);
    goalpos[0] = meta.pos_y( goal_row );
    goalpos[1] = poszy_goal[1];
    goalpos[2] = poszy_goal[0];

    if ( start_crosses==0 || goal_crosses==0 ) {
      throw std::runtime_error("AStar3DAlgo::findpath[error] start or goal point not a good 3D space point.");
    }

    // next, define the lattice
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec] * [usec/tick] = [cm/tick]
    float cm_per_wire = 0.3;
    float cm_per_row  = cm_per_tick*meta.pixel_height();
    float cm_per_col  = cm_per_pixel*meta.pixel_width();

    std::vector<float> cm_per_pixel(3);
    cm_per_pixel[0] = meta.pixel_height(); // ticks per row
    cm_per_pixel[1] = cm_per_col;
    cm_per_pixel[2] = cm_per_col;

    
    // This begins the code for the 3D linear fitter
    

    // we set t0 to be 0 on the lattice
    float tick0 = ( startpos[0] < goalpos[0] )  ? startpos[0] : goalpos[0]; 
    float y0    = ( startpos[1] < goalpos[1] )  ? startpos[1] : goalpos[1];
    float z0    = ( startpos[2] < goalpos[2] )  ? startpos[2] : goalpos[2];

    // Convert the initial tick value to an x-position value
    float x0 =  tick0 * cm_per_tick;

    // Place these points within a vector
    std::vector<float> initial_points = [x0, y0, z0];

    // To define the goal position, you have to use the same positions shown here but with the final positions flipped
    float tick_goal = ( startpos[0] < goalpos[0])   ? goalpos[0] : startpos[0];
    float y_goal    = ( startpos[1] < goalpos[1] )  ? goalpos[1] : startpos[1];
    float z_goal    = ( startpos[2] < goalpos[2] )  ? goalpos[2] : startpos[2];

    // Convert the final tick value to an x-position value
    float x_goal    = tick_goal * cm_per_tick;

    // Place these points within a vector
    std::vector<float> final_points = [x_goal, y_goal, z_goal];

    // Find the distance between the starting point and the end point 
    float distance_between_trajectory_points = sqrt((x_goal - x0)*(x_goal - x0) + (y_goal - y0)*(y_goal - y0) + (z_goal - z0)*(z_goal - z0));

    // This ends the information that can be encapsulated into the function that finds the path between the starting point and the ending point.
    Output_From_Points_On_Track     = pointsOnTrack(img_v, badch_v, initial_points, final_points, float step_size = 0.3, distance_between_trajectory_points, float min_ADC_value = 10, float neighborhood_square = 5);

    // Calculate the 'final_points_dummy' here just as a place filler; they can be in any direction from the final points from the previous function, and they should each be less than one step size away from the first set of final points so that the function is not evaluated at these final points.

    // Unfortunately I have to calculate the step sizes for each of the variables here so that I can do this.  The value of '0.5' is arbitrarily chosen, but the important point is that it is less than '1' so the full step size will not be used.
    
    std::vector<float> final_points_pseudo = [x_goal + (0.5*step_size)*(x_goal - x0)/distance_between_trajectory_points, y_goal + (0.5*step_size)*(y_goal - y0)/distance_between_trajectory_points, z_goal + (0.5*step_size)*(z_goal - z0)/distance_between_trajectory_points];

    float dist_between_end_points = sqrt(((0.5*step_size)*(x_goal - x0)/distance_between_trajectory_points) * ((0.5*step_size)*(x_goal - x0)/distance_between_trajectory_points) + ((0.5*step_size)*(y_goal - y0)/distance_between_trajectory_points) * ((0.5*step_size)*(y_goal - y0)/distance_between_trajectory_points) + ((0.5*step_size)*(z_goal - z0)/distance_between_trajectory_points) * ((0.5*step_size)*(z_goal - z0)/distance_between_trajectory_points));

    // The function 'pointsOnTrack' does not include the last point on the track's trajectory, which does not occur until after the loop breaks.  I will use the function again so that it can include the last point on the track's trajectory as well.
    Output_From_Last_Point_On_Track = pointsOnTrack(img_v, badch_v, final_points, final_points_dummy, step_size, dist_between_end_points, min_ADC_value, neighborhood_size);

    // Return the combination of 'Output_From_Points_On_Track' and 'Output_From_Last_Point_On_Track'

  }


  // Now, call a function 'pointsOnTrack' that will return all of the pixel coordinates on or near the track that should be illuminated based on how much charge the tagger finds within each of the voxels                                                                                                                                                            
  // Input values:                                                                                                                                                             
  // The image ('img_v') (type const std::vector<larcv::Image2D>&)                                                                                                              
  // The bad channel image ('badch_v') (type const std::vector<larcv::Image2D>&)                                                                                              
  // The coordinates of the initial point ('initial_point_coords') (type std::vector<float>)                                                          
  // The coordinates of the final point ('final_point_coords') (type std::vector<float>)                                                                                        
  // The step size ('step_size') (type float)                                                                                                                             
  // The total distance between the trajectory points ('distance_between_trajectory_points') (type float)                                                                 
  // The minimum ADC value ('min_ADC_value') (type float)                                                                                                                 
  // The neighborhood size that you are interested in considering ('neighborhood_square') (type 'int')
  // Now I can start writing the 'pointsOnTrack' function, which will return the number of points with charge on the track, the number of points without charge on the track, and information about neighboring points that contain charge. 

  // Note: I will consider a channel to be dead if it has a value greater than 0.0 in 'badch_v', which is the vector of plane images for the event

  std::vector<PointInfo_t> AStar3DAlgo_3DFitter::pointsOnTrack(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector<float> initial_point_coords, std::vector<float> final_point_coords, int step_size, float total_distance_between_points, float min_ADC_value, int neighborhood_size) {

    // Declare a variable equal to the 'cm_per_tick', because I need this information to convert between the x-coordinate value and the tick in the image of charge
    float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5; // [cm/usec] * [usec/tick] = [cm/tick]

    // The first step is to unpack the starting and final points
    // The 'initial_point_coords' and 'final_point_coords' are vectors with three dimensions in the detector, so I can take advantage of that in setting variables equal to these points
    float x0     =     initial_point_coords[0];
    float y0     =     initial_point_coords[1];
    float z0     =     initial_point_coords[2];

    float x_goal =     final_point_coords[0];
    float y_goal =     final_point_coords[1];
    float z_goal =     final_point_coords[2];

    // Calculate the cosine along each of the axes so that you can find the coordinates of each step that you take along the track   
    // The prescription for the cosine values is the total component of the vector in that direction (the goal point - the starting point) divided by 'distance_between_trajectory_points', the length of the particle's entire trajectory                                                                                                                     
    float cos_x = (x_goal - x0)/distance_between_trajectory_points;
    float cos_y = (y_goal - y0)/distance_between_trajectory_points;
    float cos_z = (z_goal - z0)/distance_between_trajectory_points;

    // With this information, you can calculate the step size in each dimension based on the total step size and the cosine of angle between each component and the hypotenuse of the right triangle                                                                                                                                                                 
    float x_step_size = step_size*cos_x;
    float y_step_size = step_size*cos_y;
    float z_step_size = step_size*cos_z;

    // To optimize the algorithm, I declare all of the new variables here and then fill them with each iteration of the track
    
    // Declare a vector for the numbers of the wires that correspond to a set of coordinates (only used for the y-coordinate and the z-coordinate)
    std::vector<float> wire_id(3, 0.0);
      
    // Declare the UNCOMPRESSED 'tick_value', the x-axis coordinate of the tagger before the compression factor is folded in. 
    float tick_value = 0.0;
   
    // From the step size, compute the total number of points that are located on the track.  This number will be almost certainly be a float, because the constraint on the 'step_size' as a fixed parameter does nothing to set the number of points (listed here as 'num_of_steps') to be an integer value.  I will include a test for the endpoint of the track's trajectory, which will not be included in one of the points along the track's trajectory that will undergo a test.
    float num_of_steps = total_distance_between_points/step_size;
    
    // Set 'x0', 'y0', and 'z0' equal to 'current_x_coord', 'current_y_coord', and 'current_z_coord', respectively, because these will be the variables which will be incremented
    float current_x_coord = x0;
    float current_y_coord = y0;
    float current_z_coord = z0;

    // Similarly, I will declare a variable for good pixels for each of the planes.
    int num_good_u_plane_pixels = 0;
    int num_good_v_plane_pixels = 0;
    int num_good_y_plane_pixels = 0;

    // Likewise, create a variable for the number of bad pixels found by this function.
    int num_bad_points = 0;

    // Similarly, I will declare a variable for bad pixels for each of the planes.
    int num_bad_u_plane_pixels = 0;
    int num_bad_v_plane_pixels = 0;
    int num_bad_y_plane_pixels = 0;

    // The iterator over the 'num_of_steps' can be of type 'size_t', because we are checking points less than the distance of one wire from another.
    for (size_t wire_iterator = 0; wire_iterator < num_of_steps; wire_iterator++) {

      // I will declare a data product that will be used within the 'WireCoordinate' function, just as Taritree did it
      Double_t xyz[3] {0.0, current_y_coord, current_z_coord};

      // Convert the coordinates of the wire at this value to a wire index value that can be rounded using the rounding function
      // Find the initial wire number based on the starting coordinate of the particle's path                                              
      
      // Declare a vector of the three coordinates within the detector.  I will continue down this route, because I am not sure that all of the starting points come from the same set of channels.  This allows me to do it myself instead of going back into the 'wireIntersection' function.                                                             
      // I will also round the 'wireid' value with the 'wireRoundingFunction' within the same loop.                                                      
      for (size_t wire_plane_iterator = 0; wire_plane_iterator <= 2; wire_plane_iterator++) {
	
	wire_id[wire_plane_iterator] = larutil::Geometry::GetME()->WireCoordinate( xyz , wire_plane_iterator );
	wire_id[wire_plane_iterator] = wireRoundingFunction(wire_id[wire_plane_iterator]);
	
      }

      // Once you have this information, you can convert the 'wire_id' values to the columns and row in the image based on the compression factor
      // Convert the current x0 value back to the time
      tick  = current_x_coord/cm_per_tick;

      // Declare the output data types that will be 'Double_t' and will be passed as output
      std::vector<float> output_position;
      std::vector<float> output_column;
      int output_row;                  

      // Declare a variable for the number of planes that have charge above threshold
      int num_planes_with_charge_above_threshold = 0;
      
      // Declare variables for the number of pixels on each plane along the path that have charge above threshold (a running sum that includes all of the pixels that one loops over)
      int num_u_plane_pixels_above_threshold  = 0;
      int num_v_plane_pixels_above_threshold  = 0;
      int num_y_plane_pixels_above_threshold  = 0;

      // Declare variables for the number of pixels on each plane *along the path* that are dead channels
      int num_u_plane_pixels_dead_along_path = 0;
      int num_v_plane_pixels_dead_along_path = 0;
      int num_v_plane_pixels_dead_along_path = 0;

      // I will also declare variables for 'previous above threshold row' and 'previous above threshold column' for pixels on each plane (type 'int')
      int u_previous_above_threshold_row = 0;
      int u_previous_above_threshold_col = 0;
      int v_previous_above_threshold_row = 0;
      int v_previous_above_threshold_col = 0;
      int y_previous_above_threshold_row = 0;
      int y_previous_above_threshold_col = 0;

      // I will declare variables for the 'previous dead row' and 'previous dead column' for pixels on each plane (type 'int')
      int u_previous_dead_row = 0;
      int u_previous_dead_col = 0;
      int v_previous_dead_row = 0;
      int v_previous_dead_col = 0;
      int y_previous_dead_row = 0;
      int y_previous_dead_col = 0;

      // Declare boolean variables for if a specific plane has charge above threshold
      bool u_plane_has_charge_above_threshold = false;
      bool v_plane_has_charge_above_threshold = false;
      bool y_plane_has_charge_above_threshold = false;
      
      // Initialize the variables for the row and the column
      int row = 0;
      int col = 0;
      
      // Initialize the variable for the 'pixel_score' at this coordinate value
      float pixel_score = 0.0;

      // Loop through all of the points on each of the planes to see if they contain the correct amount of charge
      // Now that I have this information, I can loop through each of the planes and see if each of the pixels have above the threshold amount of charge
      for (size_t plane_index = 0; plane_index <= 2; plane_index++) {

	// Reset the variable 'num_of_planes_with_charge_above_threshold' (This changes for each step taken)
	num_of_planes_with_charge_above_threshold = 0;

	// Declare a variable for the plane that we are currently interested in evaluating
	img = img_v.at(plane_index);
	
	// Return the row and column from the tick and the wire (the row will stay the same for the three planes, but the columns will be different)
	// The information for the compression factor is implicitly kept within there
	row = img.meta().row(tick);
	col = img.meta().col(wire_id[plane_index]);

	// Find the pixel value at this coordinate value of (row, column)
	pixel_score = img.pixel(row, col);

	// Loop through the cases for each of the planes to see which one has charge above threshold and increment the variable that
	// counts the number of planes with charge above threshold

	// U-Plane
	if (pixel_score > min_ADC_value && plane_index == 0) {
	  
	  // One more plane has charge above threshold
	  num_planes_with_charge_above_threshold += 1;

	  // The uplane does have charge above threshold on it
	  u_plane_has_charge_above_threshold = true;

	  // Check that the U plane pixel row & col are not the same as the previous one that had above the threshold amount of charge.
	  if (row != u_previous_above_threshold_row && col[0] != u_previous_above_threshold_col) {

	    // Increment 'num_u_plane_pixels_above_threshold' 
	    num_u_plane_pixels_above_threshold += 1;

	    // Set 'u_previous_above_threshold_row' and 'u_previous_above_threshold_col' equal to 'row' and 'col[0]', respectively
	    u_previous_above_threshold_row = row;
	    u_previous_above_threshold_col = col[0];
	  
	  }
	    
	}

	// V-Plane 
	if (pixel_score > min_ADC_value && plane_index == 1) {

	  // One more plane has charge above threshold
	  num_planes_with_charge_above_threshold += 1;
	  
	  // The vplane does have charge above threshold on it
	  v_plane_has_charge_above_threshold = true;

	  // Check that the V plane pixel row & col are not the same as the previous one that had above the threshold amount of charge.
	  if (row != v_previous_above_threshold_row && col[1] != v_previous_above_threshold_col) {

	    // Increment 'num_v_plane_pixels_above_threshold'
	    num_v_plane_pixels_above_threshold += 1;

	    // Set 'v_previous_above_threshold_row' and 'v_previous_above_threshold_col' equal to 'row' and 'col[1]', respectively
	    previous_above_threshold_row = row;
	    v_previous_above_threshold_col = col[1];

	  }

	}

	// Y-Plane
	if (pixel_score > min_ADC_value && plane_index == 2) {

	  // One more plane has charge above threshold
	  num_planes_with_charge_above_threshold += 1;

	  // The yplane does have charge above threshold on it
	  y_plane_has_charge_above_threshold = true;

	  // Check that the Y plane pixel row & col are not the same as the previous one that had above the threshold amount of charge.
	  if (row != y_previous_above_threshold_row && col[2] != y_previous_above_threshold_col) {

            // Increment 'num_y_plane_pixels_above_threshold'                                                                                                               
	    num_y_plane_pixels_above_threshold += 1;

            // Set 'y_previous_above_threshold_row' and 'y_previous_above_threshold_col' equal to 'row' and 'col[1]', respectively                                               
	    y_previous_above_threshold_row = row;
            y_previous_above_threshold_col = col[2];
	    
          }

	}
	
      }

	// Check to see if all of the pixels have charge above threshold.  If they do, then you do not have to search if one of the pixels is dead/look for charge in the vicinity of the main pixel
	if (num_planes_with_charge_above_threshold == 3) {
	  
	  // Add this information to the structure of type 'Double_t' in order to pass it to the next function
	  // More on all of this later
	  output_position = [current_x_coord, current_y_coord, current_z_coord];
	  output_column   = col;
	  output_row      = row;
	 
	}


	// Check to see if two of the planes have charge above threshold.  If they do, then you can diagnose the status of the channel below threshold (is it a dead channel, near dead channels, or near charge?)
	else if (num_planes_with_charge_above_threshold == 2) {

	  // Declare a variable for the 'plane_with_pixel_below_threshold'
	  int plane_with_pixel_below_threshold;

	  // Declare a variable for if the channel is a dead channel
	  bool is_lone_below_threshold_channel_dead_channel = false;

	  // If 'u_plane_has_charge_above_threshold' is false, then this is the only dead channel if this loop was entered.
	  if (u_plane_has_charge_above_threshold == false) {

	    // Set 'plane_with_pixel_below_threshold' to 0, the index for the u-plane
	    plane_with_pixel_below_threshold = 0;

	    // Find the corresponding row and column from the 'img_v.at(0)' entry, which I know has that capability.
	    row = img_v.at(0).meta().row( tick );
	    col = img_v.at(0).meta().col( wire );

	    // See if this channel is a dead channel
	    if (badch_v.at(0).pixel(row, col) > 0.0) {

	      is_lone_below_threshold_channel_dead_channel = true;

	      // Check to make sure that this dead channel was not the previous one tagged
	      if (row != u_previous_dead_row && col[0] != u_previous_dead_col) {

		// Increment 'num_u_plane_pixels_dead_along_path'
		num_u_plane_pixels_dead_along_path += 1;

		// Set 'u_previous_dead_row' and 'u_previous_dead_col' equal to the current value of the row and column
		u_previous_dead_row = row;
		u_previous_dead_col = col[0];

		// You can add the 'output_position', 'output_column' and 'output_row' to the output, the list of 'good' pixels along the path
		output_position = [current_x_coord, current_y_coord, current_z_coord];
		output_column   = col;
		output_row      = row;

	      }

	    }

	  }

	    // If 'v_plane_has_charge_above_threshold' is false, then this is the only dead channel if this loop was entered.                                         
	    if (v_plane_has_charge_above_threshold == false) {

	      // Set 'plane_with_pixel_below_threshold' to 1, the index for the v-plane                                                                                     
	      plane_with_pixel_below_threshold = 1;

	      // Find the corresponding row and column from the 'img_v.at(1)' entry, which I know has that capability.                                                     
	      row = img_v.at(1).meta().row( tick );
	      col = img_v.at(1).meta().col( wire );

	      // See if this channel is a dead channel                                                                                                                    
	      if (badch_v.at(1).pixel(row, col) > 0.0) {

		is_lone_below_threshold_channel_dead_channel = true;

		// Check to make sure that this dead channel was not the previous one tagged
		if (row != v_previous_dead_row && col[1] != v_previous_dead_col) {

		  // Increment 'num_v_plane_pixels_dead_along_path'
		  num_v_plane_pixels_dead_along_path += 1;

		  v_previous_dead_row = row;
		  v_previous_dead_col = col[1];

		  output_position = [current_x_coord, current_y_coord, current_z_coord];
		  output_column   = col;
		  output_row      = row;

		}

	      }

	    }

	      // If 'u_plane_has_charge_above_threshold' is false, then this is the only dead channel if this loop was entered.                                       
	      if (y_plane_has_charge_above_threshold == false) {

		// Set 'plane_with_pixel_below_threshold' to 2, the index for the y-plane                                                                                        
		plane_with_pixel_below_threshold = 2;

		// Find the corresponding row and column from the 'img_v.at(2)' entry, which I know has that capability.                                                  
		row = img_v.at(2).meta().row( tick);
		col = img_v.at(2).meta().col( wire );

		// See if this channel is a dead channel                                                                                                                   
		if (badch_v.at(2).pixel(row, col) > 0.0) {

		  is_lone_below_threshold_channel_dead_channel = true;

		  // Check to make sure that this dead channel was not the previous one tagged
		  if (row != y_previous_dead_row && col[2] != y_previous_dead_col) {

		    // Increment 'num_y_plane_pixels_dead_along_path'
		    num_y_plane_pixels_dead_along_path += 1;

		    y_previous_dead_row = row;
		    y_previous_dead_col = col[2];

		    output_position = [current_x_coord, current_y_coord, current_z_coord];
		    output_column   = col;
		    output_row      = row;

		  }

		}

	      }

	      // Now, if the last channel is not dead, then I should look in the vicinity of the track for a neighboring pixel with charge
	      if (is_lone_below_threshold_channel_dead_channel == false) {

		bool does_empty_pixel_have_charge_in_vicinity = doesNeighboringPixelHaveCharge(img_v, plane_with_pixel_below_threshold, row, col, neighborhood_size, min_ADC_value);
		bool does_empty_pixel_have_dead_pixel_in_vicinity = isNeighboringPixelDead(badch_v, plane_with_pixel_below_threshold, row, col, neighborhood_size);


		// If either of these evaluate to 'true', then this is a good set of pixels.  I will save the information in the output vector form.
		if (does_empty_pixel_have_charge_in_vicinity == True or does_empty_pixel_have_dead_pixel_in_vicinity == True) {

		  output_position = [current_x_coord, current_y_coord, current_z_coord];
		  output_column   = col;
		  output_row      = row;

		}

	      }
	
	} // End of test if two of the pixels have charge above threshold

	// Check to see if one of the planes have charge above threshold.  If it does, then you can diagnose the status of the two channels with charge below threshold (is it a dead channel, near dead channels, or near charge?)
	else if (num_planes_with_charge_above_threshold == 1) {

	  // Use three 'if' loops to determine which of the planes have charge below threshold
	  std::vector<int> charge_below_threshold_plane_num(2,0.0);

	  if (u_plane_has_charge_above_threshold == false && v_plane_has_charge_above_threshold == false) {

	    // Set the two values of the 'charge_below_threshold_num' equal to 0 (for the u-plane) and 1 (for the v-plane) respectively
	    charge_below_threshold_plane_num[0] = 0;
	    charge_below_threshold_plane_num[1] = 1;

	    // The 'num_y_plane_pixels_above_threshold' variable should be incremented if this is not the last 'row' and 'col' considered
	    if (row != y_previous_above_threshold_row && col[2] != y_previous_above_threshold_col) {
	      
	      // Increment 'num_y_plane_pixels_above_threshold' 
	      num_y_plane_pixels_above_threshold += 1;

	      // Set 'y_previous_above_threshold_row' and 'y_previous_above_threshold_col' equal to 'row' and 'col'
	      y_previous_above_threshold_row = row;
	      y_previous_above_threshold_col = col[2];

	    }
	      

	  }

	  if (u_plane_has_charge_above_threshold == false && y_plane_has_charge_above_threshold == false) {

	    // Set the two values of the 'charge_below_threshold_num' equal to 0 (for the u-plane) and 2 (for the y-plane) respectively                                        
	    charge_below_threshold_plane_num[0] = 0;
            charge_below_threshold_plane_num[1] = 2;

	    // The 'num_v_plane_pixels_above_threshold' variable should be incremented if this is not the last 'row' and 'col' considered                                  
	    if (row != v_previous_above_threshold_row && col[1] != v_previous_above_threshold_col) {
	      
	      // Increment 'num_y_plane_pixels_above_threshold'                                                                                                              
	      num_v_plane_pixels_above_threshold += 1;

	      // Set 'v_previous_above_threshold_row' and 'v_previous_above_threshold_col' equal to 'row' and 'col[1]'                                                       
	      v_previous_above_threshold_row = row;
	      v_previous_above_threshold_col = col;

            }

	  }

	  if (v_plane_has_charge_above_threshold == false && y_plane_has_charge_above_threshold== false) {

            // Set the two values of the 'charge_below_threshold_num' equal to 1 (for the v-plane) and 2 (for the y-plane) respectively            
	    charge_below_threshold_plane_num[0] = 1;
            charge_below_threshold_plane_num[1] = 2;

	    // The 'num_u_plane_pixels_above_threshold' variable should be incremented if this is not the last 'row' and 'col[0]' considered                            
	    if (row != u_previous_above_threshold_row && col[0] != u_previous_above_threshold_col) {

              // Increment 'num_u_plane_pixels_above_threshold'                                                                                                
	      num_u_plane_pixels_above_threshold += 1;

              // Set 'u_previous_above_threshold_row' and 'u_previous_above_threshold_col' equal to 'row' and 'col[0]'                                               
	      u_previous_above_threshold_row = row;
              u_previous_above_threshold_col = col[0];
	    }

	  }

	  // Use a vector of booleans to see if either of the pixels on the planes are dead channels
	  std::vector<bool> are_channels_dead_channels(2,false);

	  // Use a 'for' loop to check if either of the channels are dead channels
	  // The way that this works is by matching up the index for the plane number in 'charge_below_threshold_plane_num' with the index in 'are_channels_dead_channels'.
	  for (size_t dead_channel_plane = 0; dead_channel_plane <= 1; dead_channel_plane++) {

	    // Test to see if the channel is a dead channel
	    if (badch_v.at(charge_below_threshold_plane_num[dead_channel_plane]).pixel(row, col[dead_channel_plane]) > 0.0) {

	      // If the channel is dead, then set the boolean value in the 'are_channels_dead_channels' vector to true
	      are_channels_dead_channels[dead_channel_plane] = true;

	    }

	  }

	  // Declare a variable that will count how many of the channels in the planes are either (1) dead, (2) close to dead channels, or (2) close to charge if they are below threshold
	  int channels_below_threshold_dead_or_with_activity_nearby = 0;

	  // Set the boolean for if there is charge in the vicinity of a pixel to 'false'
	  bool does_empty_pixel_have_charge_in_vicinity_2D     = false;
	  
	  // Set the boolean for if there are dead pixels in the vicinity of the central pixel to 'false'
	  bool does_empty_pixel_have_dead_pixel_in_vicinity_2D = false;

	  // Of the two channels that are below threshold, if they are not dead channels, then I would like to search in their vicinity to see if there is charge (I could search in the vicinity of the dead channels as well, but I will implement that later)
	  // Start a loop over the values of 'are_channels_dead_channels', and this will lead into a loop over the 'charge_below_threshold_plane_num'
	  for (size_t dead_channel_index = 0; dead_channel_index <= 1; dead_channel_index++) {

	    if (are_channels_dead_channels.at(dead_channel_index) == true) {

	      // Check to ensure that this was not the previous pixel that was tagged as being dead
	      
	      // I will have to go through this plane-by-plane, unfortunately, based on the way that I wrote the above code.  This will be something to encapsulate in the future.
	      // U-Plane
	      if (charge_below_threshold_plane_num[dead_channel_index] == 0) {

		// Check to make sure that this is not the previous tagged dead row or dead column
		if (row != u_previous_dead_row && col[0] != u_previous_dead_col) {

		  // Increment 'num_u_plane_pixels_dead_along_path'                                                                                                            
		  num_u_plane_pixels_dead_along_path += 1;

		  // Increment the 'good_channels_above_threshold_dead_or_with_activity_nearby' 
		  good_channels_above_threshold_dead_or_with_activity_nearby += 1;

		  // Set the value of 'u_previous_dead_row' and 'u_previous_dead_col' to 'row' and 'col', respectively
		  u_previous_dead_row = row;
		  u_previous_dead_col = col;

		}

	      }

		// V-Plane
		if (charge_below_threshold_plane_num[dead_channel_index] == 1) {

		  // Check to make sure that this is not the previous tagged dead row or dead column
		  if (row != v_previous_dead_row && col[1] != v_previous_dead_col) {

		    // Increment 'num_v_plane_pixels_dead_along_path'
		    num_v_plane_pixels_dead_along_path += 1;

		    // Increment the 'good_channels_above_threshold_dead_or_with_activity_nearby' variable
		    good_channels_above_threshold_dead_or_with_activity_nearby += 1;

		    // Set the value of 'v_previous_dead_row' and 'v_previous_dead_col' to 'row' and 'col', respectively
		    v_previous_dead_row = row;
		    v_previous_dead_col = col;

		  }

		}

		// Y-Plane 
		if (charge_below_threshold_plane_num[dead_channel_index] == 2) {

		  // Check to make sure that this is not the previous tagged dead row or dead column
		  if (row != y_previous_dead_row && col[2] != y_previous_dead_col) {

		    // Increment 'num_y_plane_pixels_dead_along_path'
		    num_y_plane_pixels_dead_along_path += 1;

		    // Increment the 'good_channels_above_threshold_dead_or_with_nearby_activity' variable
		    good_channels_above_threshold_dead_or_with_activity_nearby += 1;

		    // Set the values of 'y_previous_dead_row' and 'y_previous_dead_col' to 'row' and 'col', respectively
		    y_previous_dead_row = row;
		    y_previous_dead_col = col;

		  }

		}

	    }

	    // If the channel is not a dead channel, then we need to look in its vicinity for charge/dead channels
	    else {

	      // Check to see if there is charge in the vicinity of the central pixel
	      does_empty_pixel_have_charge_in_vicinity_2D     = doesNeighboringPixelHaveCharge(img_v, charge_below_threshold_plane_num[dead_channel_index], row, col, neighborhood_size, min_ADC_value); 

	     // Check to see if there are dead pixels in the vicinity of the central pixel
	      does_empty_pixel_have_dead_pixel_in_vicinity_2D = isNeighboringPixelDeadPixel(badch_v, charge_below_threshold_plane_num[dead_channel_index], row, col, neighborhood_size);

	      // If either one of these booleans are true, then I can increment 'good_channels_above_threshold_dead_or_with_activity_nearby'
	      if (does_empty_pixel_have_charge_in_vicinity_2D == true || does_empty_pixel_have_dead_pixel_in_vicinity_2D == true) {

		good_channels_above_threshold_dead_or_with_activity_nearby += 1;

	      }

	    }

	  }

	  // If 'good_channels_above_threshold_dead_or_with_activity_nearby' is equal to 2, then this is a good set of pixels
	  if (good_channels_above_threshold_dead_or_with_activity_nearby == 2) {

	    // Add this information to the structure of type 'Double_t' in order to pass it to the next function                                                       
	    // More on all of this later                                                                                                                              
	    output_position = [current_x_coord, current_y_coord, current_z_coord];
	    output_column   = col;
	    output_row      = row;
	  
	  }

	}

	// Check to see if none of the planes have charge above threshold.  If it does, then you can diagnose of all three channels with charge below threshold (is it a dead channel, near dead channels, or near charge?)                                                                                                                                      
	else {

	  // I will declare a variable, 'good_channels_above_threshold_dead_or_with_activity_nearby_next', which mirrors the variable declared in the last loop
	  int good_channels_above_threshold_dead_or_with_activity_nearby_next = 0;

	  // Declare a vector of bools for the planes that have a dead pixel at this position
	  std::vector<bool> are_channels_dead_channels_3D(3, false);
	  
	  // Declare two booleans that will determine, if the pixel is not dead, (1) if there is charge in the vicinity and (2) if there are dead channels in the vicinity
	  // These variables will have the same name as the previous variables that served the same purpose, but they will also include the suffix '3D' for distinction

	  bool does_empty_pixel_have_charge_in_vicinity_3D     = false;
	  bool does_empty_pixel_have_dead_pixel_in_vicinity_3D = false;      

	  // First I will loop through the pixels to see if any of them are dead.
	  for (size_t dead_pixel_plane_3D = 0; dead_pixel_plane_3D <= 2; dead_pixel_plane_3D++) {

	    // Check to see if the pixel is dead
	    if (badch_v.at(dead_pixel_plane_3D).meta().pixel(row, col[dead_pixel_plane_3D]) > 0.0) {

	      // If the channel is bad, then set the entry at that index within 'are_channels_dead_channels_3D' to 'true'
	      are_channels_dead_channels[dead_pixel_plane_3D] = true;

	      // The 'if' statements for each of the 3 planes badly needs to be encapsulated.....this will be done next.
	      
	      // Include a loop over which plane it is with the dead pixels so that the number of dead pixels on each plane can be counted

	      // U-Plane
	      if (dead_pixel_plane_3D == 0) {

		// Check to make sure that this is not the same row and the same column that was last tagged as a dead pixel in the U-plane
		if (row != u_previous_dead_row && col[0] != u_previous_dead_col) {

		  // Increment the number of dead U-plane pixels found along the track's path
		  num_u_plane_pixels_dead_along_path += 1;

		  // Set 'u_previous_dead_row' to 'row' and 'u_previous_dead_col' to 'col[0]'
		  u_previous_dead_row = row;
		  u_previous_dead_col = col[0];

		}

	      }

	      // V-Plane
	      if (dead_pixel_plane_3D == 1) {

		if (row != v_previous_dead_row && col[1] != v_previous_dead_col) {

		  // Find the number of V-plane pixels found along the track's path
		  num_v_plane_pixels_dead_along_path += 1;
	      
		  // Set 'v_previous_dead_row' to 'row' and 'u_previous_dead_col' to 'col[1]'                                                                                     
		  v_previous_dead_row = row;
		  v_previous_dead_col = col[1];

		}

	      }

	      // Y-Plane
	      if (dead_pixel_plane_3D == 2) {

		if (row != y_previous_dead_row && col[2] != y_previous_dead_col) {

		  // Find the number of Y-plane pixels found along the track's path
		  num_y_plane_pixels_dead_along_path += 1;

		  // Set 'y_previous_dead_row' to 'row' and 'y_previous_dead_col' to 'col[2]'
		  y_previous_dead_row = row;
		  y_previous_dead_col = col[2];

		}
		
	      }
	      
	      // Increment 'good_channels_above_threshold_dead_or_with_activity_nearby_next', because this is a dead channel and we are not have to looking for charge in its vicinity in this version of the algorithm
	      good_channels_above_threshold_dead_or_with_activity_nearby_next += 1;
	  
	    }

	    // If the value of 'are_channels_dead_channels' at this entry was not changed to 'true', then you should check to see if there is charge or if there are dead pixels in the vicinity of the central pixel
	    if (are_channels_dead_channels[dead_pixel_plane_3D] == false) {

	      // Set 'does_empty_pixel_have_charge_in_vicinity_3D' equal to the output of the function that searches for charge in the neighborhood of the central pixel, 'doesNeighboringPixelHaveCharge', for this particular central pixel
	      does_empty_pixel_have_charge_in_vicinity_3D = doesNeighboringPixelHaveCharge(img_v, dead_pixel_plane_3D, row, col, neighborhood_size, min_ADC_value);

	      // Set 'does_empty_pixel_have_dead_pixel_in_vicinity_3D' equal to the output of the function that searches for dead pixels in the neighborhood of the central pixel, 'isNeighboringPixelDeadPixel'
	      does_empty_pixel_have_dead_pixel_in_vicinity_3D = isNeighboringPixelDeadPixel(badch_v, dead_pixel_plane_3D, row, col, neighborhood_size);

	      // If either of these booleans evaluate to 'true', then you can increment 'good_channels_above_threshold_dead_or_with_activity_nearby_next', because this is what is considered a 'good' point
	      if (does_empty_pixel_have_charge_in_vicinity_3D == true || does_empty_pixel_have_dead_pixel_in_vicinity_3D == true) {
	    
		// Increment 'good_channels_above_threshold_dead_or_with_activity_nearby_next', because this would correspond to a good pixel that we can use for our analysis
		good_channels_above_threshold_dead_or_with_activity_nearby_next += 1;

	      } // End of the loop over the 'if' statement for charge or dead pixels in the vicinity
	  
	      
	    } // End of the loop over the case in which the plane pixel being considered is not dead

	    
	  } // End of the loop over all of the planes, which in this case occurs because none of the planes have charge above threshold
      
	  // Check to see if 'good_channels_above_threshold_dead_or_with_activity_nearby_next' is equal to 3.  This means that all of the pixels are either (1) dead, (2) close to dead pixels, or (3) close to charge.  This means that we can output the position as a 'good' position.
	  if (good_channels_above_threshold_dead_or_with_activity_nearby_next == 3) {
	
	    // Add this information to the structure of type 'Double_t' in order to pass it to the next function                                                                     
	    // More on all of this later                                                                                                                                             
	    output_position = [current_x_coord, current_y_coord, current_z_coord];
	    output_column   = col;
	    output_row      = row;

	  }

	} // End the case in which none of the planes have charge above threshold ('num_planes_with_charge_above_threshold == 0')

	// Before the loop over the wires ends, I will increment 'num_of_steps'
	num_of_steps += 1;

	// Increment the coordinates being used as well so that they can be converted into wire values, which can then be converted into row values
	current_x_coord += x_step_size;
	current_y_coord += y_step_size;
	current_z_coord += z_step_size;
       
    }  // End of the loop over the wires (for some reason I'm getting mismatched parentheses)

    // I am not sure how to return the data types in a structure, but these are the quantities that I want to return
    
    // For each specific 'good point' on the tracks....
    // 'output_position'
    // 'output_column'
    // 'output_row'
    
    // For the entire track, I want to output the following quantities....
    
    // The number of pixels above threshold on each plane along the track:
    // 'num_u_plane_pixels_above_threshold'
    // 'num_v_plane_pixels_above_threshold'
    // 'num_y_plane_pixels_above_threshold'

    // The number of dead pixels on each plane along the track:
    // 'num_u_plane_pixels_dead_along_path'
    // 'num_v_plane_pixels_dead_along_path'
    // 'num_y_plane_pixels_dead_along_path'

    // Also return the final point, so that I can run the function again 
 
  } // End of the function


  // Declare a function that will search for charge in the vicinity of a pixel (one that is either dead or empty)
  // This function will be a boolean and will require the following information:
  // 'img_v' - This is the vector of images that is being examined for charge
  // 'plane_num' - This is an 'int' which corresponds to the plane that we are looking at 
  // 'central_row'       - This is an 'int' which corresponds to the row number of the central pixel
  // 'central_col'    - This is an 'int' which corresponds to the column number of the central pixel around which we will search for charge
  // 'neighborhood size' - This is an 'int' which specifies the size of the square around the central pixel that we will search for charge
  // 'min_ADC_value' - This is the minimum ADC value that a pixel in the vicinity of the central pixel can have so that it is deemed true that there is charge in the vicinity of the central pixel
  bool AStarLinearAlgo::doesNeighboringPixelHaveCharge(const std::vector<larcv::Image2D>& img_v, int plane_num, int central_row, int central_col, int neighborhood_size, float min_ADC_value) {

    // Declare a boolean, 'charge_in_vicinity', which will return 'true' if there is charge in the vicinity of the central pixel
    bool charge_in_vicinity = false;

    // You can define variables for the initial row and the final row based on where in the TPC they are
    last_row_index = img_v.front().meta().rows() - 1;
    last_col_index = img_v.front().meta().cols() - 1;

    // Declare variables for the starting and ending coordinates in the search
    int starting_row_coord = central_row - neighborhood_size;
    int ending_row_coord   = central_row + neighborhood_size;
    int starting_col_coord = central_col - neighborhood_size;
    int ending_col_coord   = central_col + neighborhood_size;

    // If I am close to the boundary of the image, then I'll have to adjust the search for this neighborhood to end on the boundary of the image.
    if (starting_row_coord < 0)
	starting_row_coord = 0;

    if (starting_row_coord > last_row_index)
      starting_row_coord = last_row_index;

    if (ending_row_coord < 0)
      ending_row_coord = 0;

    if (ending_row_coord > last_row_index)
      ending_row_coord = last_row_index;

    if (starting_col_coord < 0)
      starting_col_coord = 0;

    if (starting_col_coord > last_col_index)
      starting_col_coord = last_col_index;

    if (ending_col_coord < 0)
      ending_col_coord = 0;

    if (ending_col_coord > last_col_index)
      ending_col_coord = last_col_index;

    
    // Begin a loop over the rows and the columns of the image that are in the vicinity of 'central_row' and 'central_col' to see if there is any charge in the vicinity
    for (size_t col_iter = starting_col_coord; col_iter <= ending_col_coord; col_iter++) {

      for (size_t row_iter = starting_row_coord; row_iter <= ending_row_coord; row_iter++) {

	// First check to see if this is the same value as the central pixel.  If it is, then continue.  This is not an issue when the algorithm is looking in the vicinity of dead wires and pixels with charge below threshold, but I may change things in the future to search around pixels with charge.  Searching the central pixel is not in the spirit of the algorithm anyway.

	if (row_iter == central_row && col_iter == central_col) 
	  continue;

	// See if any pixel in this neighborhood has charge above threshold.  If it does, then return 'true' for 'charge_in_vicinity'.
	if (img_v.at(plane_num).meta().pixel(row_iter, col_iter) > min_ADC_value) {

	  charge_in_vicinity = true;

	}

      }

    }

    // Return 'charge_in_vicinity'
    return charge_in_vicinity;

  }

  // Declare a channel that will search for dead pixels in the vicinity of the central pixel.  This will take the same input parameters and use the same logic as the algorithm above, except it will take 'badch_v' as an input instead of 'img_v', and it will see if any of the entries in the vicinity of the central pixel have a value in the array that's greater than 0.0.  The 'float' input argument 'min_ADC_value' is not needed for this function.
  bool AStarLinearAlgo::isNeighboringPixelDeadPixel(const std::vector<larcv::Image2D>& badch_v, int plane_num, int central_row, int central_col, int neighborhood_size) {

    // Declare a boolean, 'charge_in_vicinity', which will return 'true' if there is charge in the vicinity of the central pixel                                               
    bool is_dead_pixel_in_vicinity = false;

    // You can define variables for the initial row and the final row based on where in the TPC they are                                                                        
    last_row_index = img_v.front().meta().rows() - 1;
    last_col_index = img_v.front().meta().cols() - 1;

    // Declare variables for the starting and ending coordinates in the search                                                                                            
    int starting_row_coord = central_row - neighborhood_size;
    int ending_row_coord   = central_row + neighborhood_size;
    int starting_col_coord = central_col - neighborhood_size;
    int ending_col_coord   = central_col + neighborhood_size;

    // If I am close to the boundary of the image, then I'll have to adjust the search for this neighborhood to end on the boundary of the image.                 
    if (starting_row_coord < 0)
      starting_row_coord = 0;

    if (starting_row_coord > last_row_index)
      starting_row_coord = last_row_index;

    if (ending_row_coord < 0)
      ending_row_coord = 0;

    if (ending_row_coord > last_row_index)
      ending_row_coord = last_row_index;

    if (starting_col_coord < 0)
      starting_col_coord = 0;

    if (starting_col_coord > last_col_index)
      starting_col_coord = last_col_index;

    if (ending_col_coord < 0)
      ending_col_coord = 0;

    if (ending_col_coord > last_col_index)
      ending_col_coord = last_col_index;

    // Begin a loop over the rows and the columns of the image that are in the vicinity of 'central_row' and 'central_col' to see if there is any charge in the vicinity     
    for (size_t col_iter = starting_col_coord; col_iter <= ending_col_coord; col_iter++) {

      for (size_t row_iter = starting_row_coord; row_iter <= ending_row_coord; row_iter++) {

        // First check to see if this is the same value as the central pixel.  If it is, then continue.  This is not an issue when the algorithm is looking in the vicinity of dead wires and pixels with charge below threshold, but I may change things in the future to search around pixels with charge.  Searching the central pixel is not in the spirit of the algorithm anyway.                                                                                                                                                                
	if (row_iter == central_row && col_iter == central_col)
	    continue;

        // See if any pixel in this neighborhood has charge above threshold.  If it does, then return 'true' for 'charge_in_vicinity'.              
	if (deadch_v.at(plane_num).meta().pixel(row_iter, col_iter) > 0.0) {

          is_dead_pixel_in_vicinity = true;

        }

      }

    }

    // Return 'charge_in_vicinity'                                                                                                                                              
    return is_dead_pixel_in_vicinity;

  }
	

  int AStarLinearAlgo::wireRoundingFunction(std::vector<float> wire_vector_value) {

    // Declare a variable that you can use for the 'integer_part' of the component
    int   integer_part; // This may have to be converted into type 'double'.  We'll see.
    float fractional_part;

    // Set the fractional part equal to the output of the 'modf' function
    fractional_part = (float)modf(wire_vector_value, &integer_part);
    // Use the ternary operator to set the value of the 'integer_part' based on if 'fractional_part' is less than or equal to 0.5
    integer_part = fractional_part >= 0.5 ? (integer_part + 1) : (integer_part + 1);

    // Return the new value for 'wire_vector_value' 
    return integer_part;
    }
  
} // This should match the end of namespace larcv......


