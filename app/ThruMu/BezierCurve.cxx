#include "BezierCurve.h"

namespace larlitecv {


	BezierCurve::BezierCurve( TMatrix& control_points ) : 
		m_cpts(control_points), kB(4,4) {
		if ( control_points.GetNcols()!=4 ) {
			throw std::runtime_error("BezierCurve must have 4 columns.");
		}
		if ( control_points.GetNrows()!=3 ) {
			throw std::runtime_error("For now BezierCurve must have 3 rows.");
		}
		setB();
	}

	BezierCurve::BezierCurve( const std::vector< std::vector<float> >& control_points ) : 
		m_cpts(control_points.front().size(), 4), kB(4,4) {
		if ( control_points.size()!=4 ) {
			throw std::runtime_error("BezierCurve must be give 4 control points.");
		}
		for (size_t i=0; i<4; i++) {
			for (size_t v=0; v<control_points.at(i).size(); v++ )
				m_cpts(v,i) = control_points.at(i).at(v);
		}
		setB();
	}



	BezierCurve::~BezierCurve() {}

	void BezierCurve::setB() {
		// define Bernstein Spline Matrix
		kB(0,0) = 1;
		kB(0,1) = -3;
		kB(0,2) = 3;
		kB(0,3) = -1;
		kB(1,0) = 0;
		kB(1,1) = 3;
		kB(1,2) = -6;
		kB(1,3) = 3;
		kB(2,0) = 0;
		kB(2,1) = 0;
		kB(2,2) = 3;
		kB(2,3) = -3;
		kB(3,0) = 0;
		kB(3,1) = 0;
		kB(3,2) = 0;
		kB(3,3) = 1;		
	}

	std::vector<float> BezierCurve::operator()(float t) {
		if (t<0 || t>1 ) {
			throw std::runtime_error( "BezierCurve::()[error]. argument must be between [0,1].");
		}
		TMatrix t_v = makeT(t);
		TMatrix sp(4,1);
		sp.Mult(kB,t_v);
		TMatrix pt(m_cpts.GetNrows(),1);
		pt.Mult(m_cpts,sp);
		std::vector<float> out( m_cpts.GetNrows(), 0.0 );
		for (int r=0; r<m_cpts.GetNrows(); r++)
			out[r] = pt(r,0);
		return out;
	}

	TMatrix BezierCurve::makeT(float t) {
		TMatrix t_v(4,1);
		t_v(0,0) = 1.0;
		t_v(1,0) = t;
		t_v(2,0) = t*t;
		t_v(3,0) = t*t*t;
		return t_v;
	}

	std::vector< std::vector<float> > BezierCurve::getCurve( int npts ) {
		float step = 1.0/float(npts);
		std::vector< std::vector<float> > curve;
		for (int istep=0; istep<=npts; istep++ ){
			float t = float(istep)*step;
			std::vector<float> pt = (*this)(t);
			curve.emplace_back( std::move(pt) );
		}
		return curve;
	}

}