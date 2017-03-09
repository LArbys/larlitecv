#ifndef __BEZIER_CURVE_H__
#define __BEZIER_CURVE_H__

#include <vector>
#include "TMatrix.h"

namespace larlitecv {

	class BezierCurve {

	public:

		BezierCurve( TMatrix& control_points );
		BezierCurve( const std::vector< std::vector<float> >& control_points );
		virtual ~BezierCurve();
		std::vector< std::vector<float> > getCurve( const int npoints );
		std::vector<float> operator()(float t);
		TMatrix makeT(float t);

		TMatrix m_cpts;
		TMatrix kB;
		void setB();

	};

}

#endif