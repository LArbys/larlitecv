#ifndef UTILFUNCTIONS_CXX
#define UTILFUNCTIONS_CXX
#include "UtilFunctions.h"

namespace larlitecv {

	TVector3 get_unit_vector(const TVector3& v1, const TVector3& v2){
		// Return unit vector from v1 to v2
		double dist = distance_between_pt(v1,v2);
		double x = (v2.x()-v1.x())/dist;
		double y = (v2.y()-v1.y())/dist;
		double z = (v2.z()-v1.z())/dist;
		TVector3 ret_vec(x,y,z);
		return ret_vec;
	}

	std::vector<int> getProjectedPixel( const TVector3& tvec,
							const larcv::ImageMeta& meta,
							const int nplanes,
							const float fracpixborder ) {
		std::vector<double> reco_vertex = {(float)tvec.X(), (float)tvec.Y(), (float)tvec.Z()};
		return getProjectedPixel( reco_vertex, meta, nplanes,fracpixborder );
	}

	std::vector<int> getProjectedPixel( const std::vector<double>& pos3d,
					    const larcv::ImageMeta& meta,
					    const int nplanes,
					    const float fracpixborder ) {
	  std::vector<int> img_coords( nplanes+1, -1 );
	  float row_border = fabs(fracpixborder)*meta.pixel_height();
	  float col_border = fabs(fracpixborder)*meta.pixel_width();

	  // tick/row
	  float tick = pos3d[0]/(::larutil::LArProperties::GetME()->DriftVelocity()*::larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3) + 3200.0;
	  if ( tick<meta.min_y() ) {
	    if ( tick>meta.min_y()-row_border )
	      // below min_y-border, out of image
	      img_coords[0] = meta.rows()-1; // note that tick axis and row indicies are in inverse order (same order in larcv2)
	    else
	      // outside of image and border
	      img_coords[0] = -1;
	  }
	  else if ( tick>meta.max_y() ) {
	    if ( tick<meta.max_y()+row_border )
	      // within upper border
	      img_coords[0] = 0;
	    else
	      // outside of image and border
	      img_coords[0] = -1;
	  }
	  else {
	    // within the image
	    img_coords[0] = meta.row( tick );
	  }

	  // Columns
	  Double_t xyz[3] = { pos3d[0], pos3d[1], pos3d[2] };
	  // there is a corner where the V plane wire number causes an error
	  if ( (pos3d[1]>-117.0 && pos3d[1]<-116.0) && pos3d[2]<2.0 ) {
	    xyz[1] = -116.0;
	  }
	  for (int p=0; p<nplanes; p++) {
	    float wire = larutil::Geometry::GetME()->WireCoordinate( xyz, p );

	    // get image coordinates
	    if ( wire<meta.min_x() ) {
	      if ( wire>meta.min_x()-col_border ) {
					// within lower border
					img_coords[p+1] = 0;
	      }
	      else
					img_coords[p+1] = -1;
	    }
	    else if ( wire>=meta.max_x() ) {
	      if ( wire<meta.max_x()+col_border ) {
					// within border
					img_coords[p+1] = meta.cols()-1;
	      }
	      else
					// outside border
					img_coords[p+1] = -1;
	    }
	    else
	      // inside image
	      img_coords[p+1] = meta.col( wire );
	  }//end of plane loop

	  // there is a corner where the V plane wire number causes an error
	  if ( pos3d[1]<-116.3 && pos3d[2]<2.0 && img_coords[1+1]==-1 ) {
	    img_coords[1+1] = 0;
	  }
	  return img_coords;
	}

	bool is_close_enough(std::vector<double> pt1, std::vector<double> compare_pt, double thresh){
		bool isclose=true;
		for (int p=0;p<3;p++){
			if ((pt1[p] > compare_pt[p]+thresh) || (pt1[p] < compare_pt[p]-thresh)) isclose=false;
		}
		return isclose;
	}

	double distance_point_plane(const TVector3& pt, double a, double b, double c, double d){
		// Calculate distance between point in TVector3 and a plane defined by:
		// ax+by+cz+d = 0
		double distance;
		double d_num = a*pt.X() + b*pt.Y() + c*pt.Z() + d;
		double d_den = std::sqrt(a*a+b*b+c*c);
		distance = (d_num/d_den);
		return distance;
	}

	std::vector<double> define_plane(const TVector3& pt1, const TVector3& pt2, const TVector3& pt3) {
		// Define a plane from 3 points. Plane equation is ax+by+cz+d=0
		// return a 4 vector of {a,b,c,d}
		std::vector<double> vec1 =  {pt2.x() - pt1.x(), pt2.y()-pt1.y(), pt2.z()-pt1.z()};
		std::vector<double> vec2 =  {pt3.x() - pt1.x(), pt3.y()-pt1.y(), pt3.z()-pt1.z()};
		// Cross product components for normal vector, use that and one of the points to get plane equation
		double cross_i = vec1[1]*vec2[2] - vec1[2]*vec2[1];
		double cross_j = vec1[2]*vec2[0] - vec1[0]*vec2[2];
		double cross_k = vec1[0]*vec2[1] - vec1[1]*vec2[0];
		double a = cross_i;
		double b = cross_j;
		double c = cross_k;
		double d = cross_i*(0-pt1.x()) + cross_j*(0-pt1.y()) + cross_k*(0-pt1.z()) ;
		std::vector<double> plane_vec = {a,b,c,d};
		return plane_vec;

	}
	TVector3 get_line_plane_intersection_pt(const TVector3& pt1,const TVector3& pt2, std::vector<double> planedef, bool check_between){
		// Take pt1 and pt2 and define the line between them, and then calculate
		// the intersection of that line with the planedef (a,b,c,d)
		// I realize I don't really have to make all these doubles, but it's clearer
		// to the reader this way.
		double a1 = pt1.X();
		double a2 = pt2.X() - pt1.X() ;
		double b1 = pt1.Y();
		double b2 = pt2.Y() - pt1.Y() ;
		double c1 = pt1.Z();
		double c2 = pt2.Z() - pt1.Z() ;
		// At this point we have the line equation:
		// P(t) =  (a1 + a2*t) X + (b1 + b2*t) Y + (b1 + b2*t) Z
		// now plut that into plane equation for x,y,z to find t, then plug
		// back into the line equation:
		double A = planedef[0];
		double B = planedef[1];
		double C = planedef[2];
		double D = planedef[3];
		double t_num = 0 - A*a1 - B*b1 - C*c1 - D;
		double t_den = A*a2 + B*b2 + C*c2;
		if (t_den == 0){
			t_den =0.01;
			std::cout << "t_den is 0...\n";
		}
		double t = t_num/t_den;
		if (check_between == true){
			// Check if t is between the t for pt1 and pt2
			double t_pt1 = (pt1.X() - a1) / a2 ;
			double t_pt2 = (pt2.X() - a1) / a2 ;
			if (a2 == 0){
				a2 =0.01;
				std::cout << "a2 is 0...\n";
			}
			if ((t_pt1 < t) && (t_pt2 < t)){
				std::cout << "Plane Not Between Points\n";
			}
			if( (t_pt1 > t) && (t_pt2 > t)){
				std::cout << "Plane Not Between Points\n";
			}
		}
		double x_inter = a1+a2*t;
		double y_inter = b1+b2*t;
		double z_inter = c1+c2*t;
		TVector3 intersection_pt(x_inter,y_inter,z_inter);
		return intersection_pt;
	}

	double distance_between_pt(const std::vector<double>& pt1,const std::vector<double>& pt2){
		/*
		This function finds the distance between pt1 and pt2
		which should both be in x,y,z
		*/

	  double distance = -1;
	  distance = std::sqrt(std::pow(pt1[0]-pt2[0],2)+std::pow(pt1[1]-pt2[1],2)+std::pow(pt1[2]-pt2[2],2));
	  return distance ;
	}
	double distance_between_pt(const TVector3& pt1,const TVector3& pt2){
		/*
		This function finds the distance between pt1 and pt2
		which should both be TVector3
		*/

		double distance = -1;
		distance = std::sqrt(std::pow(pt1.x()-pt2.x(),2)+std::pow(pt1.y()-pt2.y(),2)+std::pow(pt1.z()-pt2.z(),2));
		return distance;
	}
	double distance_between_pt2d(const int x1, const int y1, const int x2, const int y2){
		/*
		This function finds the distance between pt1 and pt2
		which should both be in x,y
		*/
		double xd1 = (double)x1;
		double yd1 = (double)y1;
		double xd2 = (double)x2;
		double yd2 = (double)y2;
		double distance = -1;
		distance = std::sqrt(std::pow(xd1-xd2,2)+std::pow(yd1-yd2,2));
		return distance ;
	}

	TVector3 perform_sce_fwd(larlite::mcstep const& mcstep){
		/*
		This function takes an roi, and applies an sce_fwd correction to it,
		in order to match Reconstruction. It requires global variables g_driftv
		g_samplingrate, and g_sce_class
		-Josh
		*/
		larutil::SpaceChargeMicroBooNE g_sce_class; // default is fwd mcc9
		double g_driftv = larutil::LArProperties::GetME()->DriftVelocity();
		double g_samplingrate = larutil::DetectorProperties::GetME()->SamplingRate();

	  TVector3 xyz;
	  double tick = mcstep.X() / g_driftv * g_samplingrate / 1000 + 3200;
	  double delta_x_tickhack =  ( g_driftv ) * 1000 / g_samplingrate * ((6+tick-3200)*0.014);
		std::vector<double> xyzchange(3,0);
		if ( error_check(mcstep.X(),mcstep.Y(),mcstep.Z()) ){
			xyzchange = g_sce_class.GetPosOffsets(mcstep.X(),mcstep.Y(),mcstep.Z());
		}
	  double new_x = mcstep_time_adjust(mcstep);
	  xyz.SetX( new_x - xyzchange[0] + delta_x_tickhack + 0.6);
	  xyz.SetY( mcstep.Y() + xyzchange[1] ) ;
	  xyz.SetZ( mcstep.Z() + xyzchange[2] ) ;
	  return xyz;
	}
	double mcstep_time_adjust(const larlite::mcstep& mcstep){
		/*
		This function takes an mcstep and returns the new adjusted x position
		in accordance with Taritree's Hack
		-Josh
		*/
	  double cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5;
	  double time = mcstep.T();
	  double tick = larutil::TimeService::GetME()->TPCG4Time2Tick(time) + mcstep.X()/(cm_per_tick);
	  double new_x = (tick - 3200)*cm_per_tick;
	  return new_x;
	}
	bool error_check(double x, double y, double z){
		/*
		This file checks if x,y,z are outside the detector's dimensions as needed
		for the GetPosOffsets function, to avoid throughing warnings endlessly
		-Josh
		*/
		if ((x == 0.) || (x > 256.)){
			return false;
		}
		else if ((y <= -116.5) || (y > 116.5)){
			return false;
		}
		else if ((z == 0.) || (z > 1037.)){
			return false;
		}
		return true;
	}
}
#endif
