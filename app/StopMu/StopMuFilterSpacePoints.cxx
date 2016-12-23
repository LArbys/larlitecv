#include "StopMuFilterSpacePoints.h"


namespace larlitecv {

  StopMuFilterSpacePoints::StopMuFilterSpacePoints( const StopMuFilterSpacePointsConfig& config  ) {
    m_config = &config;
  }

  std::vector<larcv::Pixel2D> StopMuFilterSpacePoints::filterSpacePoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
									  const std::vector<larcv::Image2D>& thrumu_pixels,
									  const std::vector<larcv::Image2D>& badch_imgs) {
    // this function selects spacepoints we want to pass to the stop-mu tracker
    // the goal is to
    // (1) filter out space points that below any through-going muon points
    // (2) merge points that are on the same track (choosing the one closest to the boundary)

    std::cout << "StopMuFilterSpacePoints::StopMuFilterSpacePoints" << std::endl;

    std::vector< std::vector<const larcv::Pixel2D*> > not_thrumu;
    removeThroughGoingEndPoints( spacepoints_list, thrumu_pixels, not_thrumu );

    std::cout << "number of non-thru-mu endpoints: " << not_thrumu.size() << std::endl;

    

  }
  
  void StopMuFilterSpacePoints::removeThroughGoingEndPoints( const std::vector< larcv::EventPixel2D*  >& spacepoints_list,
							     const std::vector<larcv::Image2D>& thrumu_pixels,
							     std::vector< std::vector<const larcv::Pixel2D*> >& passlist ) {
    for ( auto &ev_ptr : spacepoints_list ) {
      size_t nendpts = (*ev_ptr).Pixel2DArray(0).size();
      for (size_t i=0; i<nendpts; i++) {
	// collect the pixels across the planes
	std::vector<const larcv::Pixel2D*> pix_v(3,0);
	for (int p=0; p<3; p++) {
	  pix_v[p] = &(ev_ptr->Pixel2DArray(p).at(i));
	}
	
	// we check the thrumu_pixel, is there a thru-mu tag nearby?
	int nplanes_tagged = 0;
	for (int p=0; p<3; p++) {
	  int pix_row=pix_v[p]->Y();
	  int pix_col=pix_v[p]->X();
	  for (int r=-m_config->row_tag_neighborhood; r<=m_config->row_tag_neighborhood; r++) {
	    for (int c=-m_config->col_tag_neighborhood; c<=m_config->col_tag_neighborhood; c++) {
	      int row = pix_row+r;
	      int col = pix_col+c;
	      if ( row<0 || row>=thrumu_pixels.at(p).meta().rows() || col<0 || col>=thrumu_pixels.at(p).meta().cols() )
		continue;
	      if ( thrumu_pixels.at(p).pixel(row,col)>0 )
		nplanes_tagged++;
	    }
	  }//end of r loop
	}//end of plane loop

	if ( nplanes_tagged>=2 )
	  passlist.push_back( pix_v );

      }//end of endpoint loop
    }//end of list loop
  }//end of removeThroughGoingEndPoints
							     

}
