import os,sys
import pyqtgraph as pg
from pyqtgraph import PlotDataItem
import numpy as np
from larcv import larcv

class VisFlashHits:

    COLORS = [ (76,0,153,125), # purple
               (102,255,102),  # green
               (51,255,255),   # cyan
               (255,128,0),    # orange
               (0,76,153),     # blue
               (204,0,0),      # red
               (0,102,0),      # dark-green
               (102,0,51),     # dark-magenta
               (255,255,0),    # yellow
               ]
    NCOLORS = len(COLORS)
    PLANE_COLORS = {0:(255,204,204,125),
                    1:(51,255,51,125),
                    2:(102,102,255,125)}
    CH_MARKER = { 0:"t",
                  1:"o",
                  2:"s",
                  3:"d" }
    
    def __init__(self):
        pass

    def configure(self,pset):
        self.imgproducer = pset.get("ImgProducer")
        self.producer    = pset.get("FlashHitProducer")

    def visualize( self, larlite_io, larcv_io, rawdigit_io ):
        event_imgs   = larcv_io.get_data( larcv.kProductImage2D, self.imgproducer )
        event_pixels = larcv_io.get_data( larcv.kProductPixel2D, self.producer )

        lcv_imgs = event_imgs.Image2DArray()

        cluster_vecs = [] # output container

        i = 0
        npts = 0
        for iimg in xrange(0,lcv_imgs.size()):
            lcv_img = lcv_imgs.at(iimg)
            meta = lcv_img.meta()
            plane = meta.plane()
            pixel_v  = event_pixels.Pixel2DArray(plane)
            nends = pixel_v.size()
            x = np.zeros( nends )
            y = np.zeros( nends )
            for ihit in range(0,nends):
                pixel = pixel_v.at( ihit )
                x[ihit] = meta.pos_x(pixel.X())
                y[ihit] = meta.pos_y(pixel.Y())
                npts+=1
            color = VisFlashHits.PLANE_COLORS[ plane ]
            plot = PlotDataItem( x=x, y=y, 
                                 pen=None,
                                 symbolBrush=pg.mkBrush(color=color),
                                 symbol='+',symbolPen=pg.mkPen(color=color,width=0.0) )
            cluster_vecs.append(plot)

        print "flash plots: ",len(cluster_vecs)," number of points: ",npts
        return cluster_vecs
            
