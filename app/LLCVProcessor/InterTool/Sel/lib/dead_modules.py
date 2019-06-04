import ROOT
from larcv import larcv as lcv

def BlankImage():
    name = "BlankImage"
    return lcv.BlankImage(ROOT.std.string(name))

def WireMask():
    name = "WireMask"
    return lcv.WireMask(ROOT.std.string(name))

def attach_dead(proc):
    
    blank_image = BlankImage()
    proc.add_lc_proc(blank_image)

    wire_mask = WireMask()
    proc.add_lc_proc(wire_mask)
    
