#include "PushBoundarySpacePoint.h"

#include "UBWireTool/UBWireTool.h"
#include "TaggerTypes/dwall.h"

namespace larlitecv {

    PushBoundarySpacePoint::PushBoundarySpacePoint()
    : m_foxalgo_cfg()
    , m_foxalgo( m_foxalgo_cfg )
    {
        m_foxalgo_cfg.step_size = 4.0;
        m_foxalgo_cfg.num_step_attempts = 2;
        m_foxalgo_cfg.pixel_thresholds.resize(3,10.0);
        m_foxalgo_cfg.min_hit_width = 1;
        m_foxalgo_cfg.hit_neighborhood = 2;
        m_foxalgo_cfg.segment_frac_w_charge = 0.5;
        m_foxalgo_cfg.radius_reduction_factor = 0.5;
        m_foxalgo_cfg.min_cosine = 0.0;
        m_foxalgo_cfg.max_steps = 100;
        m_foxalgo_cfg.verbosity = 0;
        FoxTrotTrackerAlgo reconfiged(m_foxalgo_cfg);
        std::swap(m_foxalgo,reconfiged);
    }

    BoundarySpacePoint PushBoundarySpacePoint::pushPoint( const BoundarySpacePoint& sp, const std::vector<larcv::Image2D>& img_v,
        const std::vector<larcv::Image2D>& badch_v ) {
        clear();


        //std::cout << "Given point starts at : (" << sp.pos()[0] << "," << sp.pos()[1] << "," << sp.pos()[2] << ")" << std::endl;

        FoxTrack track = runFoxTrot( sp, img_v, badch_v );

        //std::cout << "Initial fox trot track ends at: (" << track.back().pos()[0] << "," << track.back().pos()[1]
        //            << "," << track.back().pos()[2] << ")" << std::endl;

	// Print 
	
        BoundarySpacePoint endpt = scanTrackForEndPoint( sp, track, img_v, badch_v );
        BoundarySpacePoint endpt_out = evalEndPoint( endpt );
        // save data for plotting or debugging
        m_tracklist_v.emplace_back( std::move(track) );
        m_endpoints_v.emplace_back( std::move(endpt) );
        m_endpoints_v.push_back( endpt_out );
        return endpt_out;
    }

    void PushBoundarySpacePoint::clear() {
        m_tracklist_v.clear();
        m_endpoints_v.clear();
    }

    FoxTrack PushBoundarySpacePoint::runFoxTrot( const larlitecv::BoundarySpacePoint& sp, const std::vector<larcv::Image2D>& img_v,
                    const std::vector<larcv::Image2D>& badch_v ) {

        std::vector<float> startdir(3,0.0);
        if ( sp.type()<larlitecv::kImageEnd ) {
            switch ( sp.type() ) {
            case 0: // top
                startdir[1] = 1.0;
                break;
            case 1: // bottom
                startdir[1] = -1.0;
                break;
            case 2: //upstream
                startdir[2] = -1.0;
                break;
            case 3: // downstream
                startdir[2] = 1.0;
                break;
            case 4: // anode
                startdir[0] = -1.0;
                break;
            case 5: // cathode
                startdir[0] = 1.0;
                break;
            default:
                throw std::runtime_error( "PushBoundarySpacePoint::runFoxTrot -- unrecognized boundary point type" );
                break;
            };
        }//end of if ! imagends
        else {
            std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( sp.pos(), img_v.front().meta(), (int)img_v.size() );
            startdir[0] = ( imgcoords[0] < (int)img_v.front().meta().rows()/2 ) ? 1.0 : -1.0;
        }
	   //std::cout << "starting foxtrot for type=" << sp.type() << " and dir=(" << startdir[0] << "," << startdir[1] << "," << startdir[2] << ")" << std::endl;
        return m_foxalgo.followTrack( img_v, badch_v, badch_v, sp, startdir );
    }

    BoundarySpacePoint PushBoundarySpacePoint::scanTrackForEndPoint( const BoundarySpacePoint& original, const FoxTrack& track,
        const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
        // we go through small steps along the path, finding the end point
        // we find the end point by taking sub steps between the fox trot steps
        // we record the last point with good charge

        if (track.size()==1) {
            return BoundarySpacePoint(original);
        }

        std::vector<float> endpos(3,0);
        std::vector<float> enddir(3,0);
        float max_step_size = 0.3;
        int hit_neighborhood = 2;
	bool pixel_in_image  = true;
	
        std::vector<float> thresholds(3,10.0);
        for ( int istep=1; istep<(int)track.size(); istep++ ) {
	  // I will make this a constant below after I make sure it is inside the image.
            const FoxStep& last_step    = track[istep-1];
            const FoxStep& current_step = track[istep];
            float dir[3] = {0.0};
            float dist = 0.0;
            for (int i=0; i<3; i++) {
                dir[i] = current_step.pos()[i]-last_step.pos()[i];
                dist += dir[i]*dir[i];
            }
            dist = sqrt(dist);
            for (int i=0; i<3; i++) dir[i] /= dist;
            int nsubsteps  = dist/max_step_size+1;
            float stepsize = dist/nsubsteps;
            int substeps_w_charge = 0;

            std::vector<float> last_good_point(3,0);
            for (int i=0; i<3; i++)
	      last_good_point[i] = last_step.pos()[i];
	      
		   
            for (int isubstep=0; isubstep<=nsubsteps; isubstep++) {
                std::vector<float> subpos(3,0.0);
                for (int i=0; i<3; i++)
                    subpos[i] = last_step.pos()[i] + isubstep*dir[i]*stepsize;
                std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( subpos, img_v.front().meta(), (int)img_v.size() );

		// Make sure that the central pixel is located in the image.
		// Row Pixel.
		if  (imgcoords[0] < 0 || imgcoords[0] >= (int)img_v.front().meta().rows())
		  pixel_in_image = false;

		// Column Pixel.
		for (size_t in_img_iter = 0; in_img_iter<3; in_img_iter++){
		  if (imgcoords[in_img_iter+1] < 0 || imgcoords[in_img_iter+1] >= (int)img_v.front().meta().cols()){
		    pixel_in_image = false; }
		}

		// Continue if 'pixel_in_image' is false - this central pixel is out of the image and will create problems.
		if (pixel_in_image == false) {
		  continue;
		}

                bool hascharge = false;
                bool found_substep_end = false;
                int ngoodplanes = 0;
                for (size_t p=0; p<img_v.size(); p++) {
		            bool foundhit_on_plane = false;
		            for (int dr=-hit_neighborhood; dr<=hit_neighborhood; dr++) {
                        int row = imgcoords[0]+dr;
                        if ( row<0 || row>=(int)img_v.front().meta().rows()) continue;
                        for (int dc=-hit_neighborhood; dc<=hit_neighborhood; dc++) {
                            int col = imgcoords[p+1] + dc;
                            if ( col<0 || col>=(int)img_v.front().meta().cols()) continue;
                            if ( img_v[p].pixel(row,col)>thresholds[p] || badch_v[p].pixel(row,col)>0 ) {
                                foundhit_on_plane = true;
                            }
                            if ( foundhit_on_plane )
                                break;
                        }// col neighborhoood
                        if ( foundhit_on_plane ) {
                            break;
                        }
		            }// end of row loop
		            if ( foundhit_on_plane )
		              ngoodplanes++;
                }//end of plane loop
                if ( ngoodplanes>=3 ) {
                    hascharge = true;
                }

                if ( hascharge ) {
                    substeps_w_charge += 1;
                    if (!found_substep_end) {
                        for (int i=0; i<3; i++)
                            last_good_point[i] = subpos[i];
                        found_substep_end = true;
                    }
                }

                // std::cout << " checking substep (" << istep << "," << isubstep << "/" << nsubsteps << "): "
                //             << "hascharge=" << hascharge << " "
                //             << "found_substep_end=" << found_substep_end
                //             << std::endl;
            }// end of substep loop

            // fill the end point with the last good subpos end
            for (int i=0; i<3; i++) {
                endpos[i] = last_good_point[i];
                enddir[i] = dir[i];
            }
            if ( float(substeps_w_charge)/(nsubsteps) < 0.9 ) {
                // this step is probably the last good step. stop here and fill the end pos
                break;
            }
        } // end of step loop

        std::vector<BoundaryEndPt> endpt_v;
        std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( endpos, img_v.front().meta(), (int)img_v.size() );
        for (int p=0; p<3; p++) {
            BoundaryEndPt endpt_t( imgcoords[0], imgcoords[p+1], original.type() );
            endpt_v.emplace_back( std::move(endpt_t) );
        }
        BoundarySpacePoint pushed( original.type(), std::move(endpt_v), endpos[0], endpos[1], endpos[2] );

        return pushed;
    }//end of find end point

    larlitecv::BoundarySpacePoint PushBoundarySpacePoint::evalEndPoint( const larlitecv::BoundarySpacePoint& moved_pt ) {
        // with the moved end point, do we need to update the type?
        int boundary_type = moved_pt.type();
        float dist = larlitecv::dwall( moved_pt.pos(), boundary_type );
        if ( (int)moved_pt.type()==boundary_type ) {
            std::vector<BoundaryEndPt> endpt_v;
            for (int p=0; p<3; p++)
                endpt_v.push_back( moved_pt[p] );
            return BoundarySpacePoint( (larlitecv::BoundaryEnd_t)boundary_type, std::move(endpt_v),
                                        moved_pt.pos()[0], moved_pt.pos()[1], moved_pt.pos()[2] );
        }
        return BoundarySpacePoint( moved_pt );
    }

}
