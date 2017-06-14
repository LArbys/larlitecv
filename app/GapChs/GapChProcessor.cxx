#include "GapChProcessor.h"
#include "EmptyChannelAlgo.h"

namespace larlitecv {

	void GapChProcessor::addEventImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
		// we make copies!
		m_evimg_list.push_back( img_v );
		m_evbadch_list.push_back( badch_v );
	}

	std::vector<larcv::Image2D> GapChProcessor::makeGapChImage() {
		std::vector<larcv::Image2D> gapch_v;

		if ( m_evimg_list.size()==0 )
			return gapch_v;

		const size_t nplanes = m_evimg_list.front().size();
		for (size_t p=0; p<nplanes; p++) {
			larcv::Image2D img( m_evimg_list.front().at(p).meta() );
			img.paint(255.0); // all channels tagged as missing, unless seen otherwise
			gapch_v.emplace_back( std::move(img) );
		}


		EmptyChannelAlgo emptyalgo;
		std::vector<int> ngood(3,0);
		for ( size_t ev=0; ev<m_evimg_list.size(); ev++ ) {
			// for each image, we generate an image showing the empty channels
			// we combine with this the cumulative gapch image
			// basically, if a channel is marked empty unless it has been seen to not be empty
			const std::vector< larcv::Image2D >& img_v   = m_evimg_list.at(ev);
			const std::vector< larcv::Image2D >& badch_v = m_evbadch_list.at(ev);

			// use the empty algo to find the missing channels
			std::cout << "nimgs in " << ev << ": " << img_v.size() << std::endl;
			std::vector<larcv::Image2D> evmissingchs = emptyalgo.findMissingBadChs( img_v, badch_v, 15, 3456, 100 );

			// compare list with existing gapchs
			for (size_t p=0; p<nplanes; p++) {
				larcv::Image2D& gapch = gapch_v.at(p);
				const larcv::ImageMeta& meta = gapch.meta();
				int nmissing = 0;
				for ( size_t c=0; c<meta.cols(); c++ ) {
					float gapchval = gapch.pixel(10,c);
					if (gapchval==0) 
						continue; // once untagged as gapch, cannot switch back
					float missingchval = evmissingchs.at(p).pixel(10,c);
					if ( missingchval>0 ) {
						nmissing++;
					}
					else {
						gapch.paint_col(c,0);
						ngood[p]++;
					}
				}
				std::cout << " plane" << p << " missing=" << nmissing << std::endl;
			}//end of loop over planes
			std::cout << "after event " << ev << ": plane0=" << ngood[0] << " plane1=" << ngood[1] << " plane2=" << ngood[2] << std::endl;
		}//end of loop over all imgs in list

		std::cout << "Number of GapChs found after incorporating " << m_evimg_list.size() << " event image sets." << std::endl;
		for (size_t p=0; p<nplanes; p++) {
			int ngapchs = 0;
			for (size_t c=0; c<gapch_v.at(p).meta().cols(); c++) {
				if ( gapch_v.at(p).pixel(0,c)>0 )
					ngapchs++;
			}
			std::cout << "  plane " << p << ": " << ngapchs << std::endl;
		}

		return gapch_v;

	}


}