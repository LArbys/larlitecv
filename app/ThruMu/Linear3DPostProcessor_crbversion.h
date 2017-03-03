#ifndef __LINEAR3DPOSTPROCESSOR_H__
#define __LINEAR3DPOSTPROCESSOR_H__

/*! \brief Post-process BMTrackCluster3D objects from Linear3DFitter
 *
 *
 * This class processes a vector of BMTrackCluster3D objects.
 * It tries to find tracks that are a subsegment of a larger track.
 * Also, if it can merge two tracks, it will try to do that too.
 *
 * Basically, it tries to clean up the output and make up for failures of 
 * the Linear3DFitter.
 *
 * author(s): Taritree Wongjirad (taritree@mit.edu)
 *
 * revisions
 * 2/14/2017: first writing
 */

#include <vector>

// larlite/UserDev/BasicTool/GeoAlgo
#include "GeoAlgo.h"
#include "GeoLineSegment.h"

// larlitecv/app/Thrumu
#include "BMTrackCluster3D.h"

namespace larlitecv {

  class Linear3DPostProcessor {

  public:
    Linear3DPostProcessor() {};
    virtual ~Linear3DPostProcessor() {};

    std::vector < BMTrackCluster3D > process( std::vector< BMTrackCluster3D >& tracks_v );

    bool pointProximity(geoalgo::Line long_line, geoalgo::LineSegment* short_linesegment, int resolution, double maximum_dist);

    int longestDimFinder(const std::vector <double> point_one_coords, const std::vector <double> point_two_coords);

    std::vector < std::vector<double> > getTrackPoints(const std::vector<double> point_one_coords, const std::vector<double> point_two_coords, int resolution, std::vector <int> shorter_dim_idx_vector, int longest_dim_idx);

    std::vector <std::vector <double> > longerTrackExtensionCoordinates(std::vector<double> longer_line_first_point_coords, std::vector<double> longer_line_second_point_coords, std::vector<double> shorter_linesegment_first_point_coords, std::vector<double> shorter_linesegment_second_point_coords);

    void  concatenatingTwoTrackParts(std::vector< BMTrackCluster3D >& tracks_v, int itrack_a, int itrack_b, BMTrackCluster3D& track_a, BMTrackCluster3D& track_b, std::vector < std::vector <double> > endpoint_coordinates_of_new_track);


  };

}

#endif
