#ifndef __COMBINATOR_H__
#define __COMBINATOR_H__

#include <vector>
#include <iostream>

/*
	Often we need to evaluate combinations of objects across several planes.
	This class takes in a vector of vectors and provides a way to iterator over combinations among elemnts of the vectors.

*/
namespace larlitecv {

	template <class T> class Combinator {
	  public:
	  	Combinator( const std::vector< std::vector<T> >& data )
	  	 : m_data(data) {
	  	 	m_combo.resize( m_data.size(), 0 );
	  	};
	  	virtual ~Combinator() {};

	  	std::vector< const T* > getCombo() {
	  		std::vector< const T* > out;
	  		if ( isLast() ) return out;
	  		for ( int p=0; p<ndims(); p++ )
  	  		out.push_back( &(m_data.at(p).at( m_combo[p] )) );
  	  	return out;
	  	}
	  	std::vector<int> getIndexCombo() {
	  		return m_combo;
	  	};
	  	bool isLast() {
	  		if ( m_combo[ndims()-1]!=dimsize(ndims()-1) )
					return false;
				return true;	  		
	  	};
	  	bool next() {
				if ( isLast() ) return true; // prevent moving further
				for ( int p=0; p<ndims(); p++) {
					m_combo[p]++;
					if ( m_combo[p]<dimsize(p) )
						break;
					else if ( p!=(ndims()-1) ) {
						m_combo[p] = 0;
					}
					// reset this value, move to next dim
				}
				return false;
	  	};

	  protected:
	  	std::vector< int > m_combo;
	  	const std::vector< std::vector<T> >& m_data;
	  	int ndims() { return m_data.size(); };
	  	int dimsize( int i ) { return m_data.at(i).size(); };
  };

}


#endif
