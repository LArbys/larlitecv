#ifndef __CLUSTERGROUPMATCHINGTYPES_H__
#define __CLUSTERGROUPMATCHINGTYPES_H__

#include <vector>
#include "ClusterGroupAlgo.h"

namespace larlitecv {

	// The object returned by ClusterGroupMatchingAlgo is ChargeVolume, which you can find later on in the file.

	// ==================================================================
	// Below is a lot of internal Types used by ClusterGroupMatchingAlgo.
	// For the most part, not needed by user.
	// Actually, if I remove the Slices_t definition, I can decouple...

	// triple for tracking indices of matched cluster group in some container
	class GroupIndices_t : public std::vector<int> {
	public:
		GroupIndices_t(int idx1, int idx2, int idx3 ) {
			resize(3,0);
			at(0) = idx1;
			at(1) = idx2;
			at(2) = idx3;
		};
		virtual ~GroupIndices_t() {};
	};


	// We use this to metric to rank preliminary matches.  
	// We use the time diff of start and end point to set the order of processing candidate cluster group matches.
	class PreMatchMetric_t {
	public:
		PreMatchMetric_t( const ClusterGroup& p1, const ClusterGroup& p2, const ClusterGroup& p3 );
		virtual ~PreMatchMetric_t() {};

		bool operator()(const PreMatchMetric_t& lhs, const PreMatchMetric_t& rhs ) { //< sorting comparator
			if ( lhs.dtSpan + lhs.dtEnd < rhs.dtSpan+rhs.dtEnd ) 
				return true;
			return false;
			// if ( lhs.dtSpan<rhs.dtSpan)
			// 	return true;
			// else if ( lhs.dtSpan==rhs.dtSpan && lhs.dtEnd<rhs.dtEnd )
			// 	return true;
			// return false;
		};

		bool operator<( const PreMatchMetric_t& rhs ) { //< sorting comparator
			return (*this)( *this, rhs );
		};

		float dtSpan; // difference in time-spane
		float dtEnd;  // difference in end-points
		//float dQ;   // difference in charge (in future?)			
		std::vector< const ClusterGroup* > m_clusters;
		std::vector< int > m_index_combo;
	};

	// Intersection Points in the (y,z) coordinates. Basically the projection on the wire planes.
	class Point_t : public std::vector<float>  {
	public:
		Point_t( float y, float z, float start, float end )
		 : tickstart(start), tickend(end) { 
		  resize(2); 
		  at(0)=y; 
		  at(1)=z; 
		};
		virtual ~Point_t() {};
		float tickstart;
		float tickend;
		const std::vector<float> as_vec() const {
			std::vector<float> vec(2);
			for (int i=0; i<2; i++) 
				vec[i] = (*this)[i];
			return vec;
		};
	};

	// Some containers for our Point_t object
  typedef std::vector<Point_t> PointList_t;
  typedef std::vector<PointList_t> Slices_t;

  // Defines a wire range: just two ints
	class WireInterval : public std::vector<int> {
	public:
		WireInterval() {
			resize(2);
			at(0) = -1;
			at(1) = -1;
		};
		WireInterval(int w1, int w2) {
			resize(2);
			at(0) = w1;
			at(1) = w2;
		};
		virtual ~WireInterval() {};
	};

	// ================================================================
	// This is the output object of ClusterGroupMatchingAlgo
	// ================================================================
	
  class ChargeVolume {
  public:
  	ChargeVolume() {
  		frac_good_slices = 0.;
  		num_slices = 0;
  		num_good_slices = 0;
  		plane_charge.resize(3,0.0);
  		clustergroup_indices.resize(3,-1);
  	};
  	virtual ~ChargeVolume() {};
  	Slices_t slices;
  	float frac_good_slices;
  	int num_good_slices;
  	int num_slices;
  	std::vector<float> plane_charge;
  	std::vector<int> clustergroup_indices;

		bool operator<(const ChargeVolume& rhs ) {
			if ( frac_good_slices>rhs.frac_good_slices)
				return true;
			return false;
		};

		bool isempty() {
			if ( clustergroup_indices[0]==-1 && clustergroup_indices[1]==-1 && clustergroup_indices[2]==-1 )
				return true;
			return false;
		};
	};	

}

#endif