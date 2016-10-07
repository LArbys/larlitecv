import os,sys
import pyqtgraph as pg
from pyqtgraph import PlotDataItem
import numpy as np
from larcv import larcv

class VisBoundaryHits:

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
    PLANE_COLORS = {0:(125,0,0,125),
                    1:(0,125,0,125),
                    2:(0,0,125,125)}
    CH_MARKER = { 0:"t",
                  1:"o",
                  2:"s",
                  3:"d" }               
    
    def __init__(self):
        pass

    def configure(self,pset):
        self.img_producer = pset.get("ImageProducer")
        self.top_producer = pset.get("top")
        self.bottom_producer = pset.get("bottom")
        self.upstream_producer = pset.get("upstream")
        self.downstream_producer = pset.get("downstream")
        self.plane_selection = pset.get("plane",-1)
        self.producers = [ self.top_producer, self.bottom_producer, self.upstream_producer, self.downstream_producer ]

    def visualize( self, larlite_io, larcv_io, rawdigit_io ):
        endpt_plots = []
        event_imgs = larcv_io.get_data( larcv.kProductImage2D, self.img_producer )
        if event_imgs.Image2DArray().size()==0:
            print "no images? stop boundary hit vis."
            return
        meta = event_imgs.Image2DArray().at(0).meta()
        
        for iproducer,producer in enumerate( self.producers ):
            
            event_pixels = larcv_io.get_data( larcv.kProductPixel2D, producer )

            for iplane in range(0,3):
                if self.plane_selection>=0 and iplane!=self.plane_selection:
                    continue

                endpts = event_pixels.Pixel2DArray( iplane )
                nendpts = endpts.size()
                if nendpts==0:
                    continue
                x = np.zeros( nendpts )
                y = np.zeros( nendpts )
                for ipt,endpt in enumerate(endpts):
                    x[ipt] = meta.pos_x( endpt.X() )
                    y[ipt] = meta.pos_y( endpt.Y() )
                color  = VisBoundaryHits.PLANE_COLORS[iplane]
                symbol = VisBoundaryHits.CH_MARKER[iproducer]
                plot = PlotDataItem( x=x, y=y, 
                                     pen=None,
                                     symbolBrush=pg.mkBrush(color=color),
                                     symbol=symbol,symbolPen=pg.mkPen(color=color,width=0.0) )
                endpt_plots.append( plot )

        print "number of end point collections: ",len(endpts)
        return endpt_plots
            
