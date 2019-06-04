import ROOT
from larcv import larcv as lcv

def ChannelMax(plane):
    name = "ChannelMaxPlane%d" % plane
    return lcv.ChannelMax(ROOT.std.string(name))

def Combine():
    name = "CombineImages"
    return lcv.CombineImages(ROOT.std.string(name))

def Mask(prefix=""):
    name = "%sSegment" % prefix
    return lcv.SegmentMask(ROOT.std.string(name))

def Image(prefix=""):
    name = "%sImage" % prefix
    return lcv.MaskImage(ROOT.std.string(name))


def attach_ssnet(proc):
    # max
    for plane in xrange(3):
        res = ChannelMax(plane)
        proc.add_lc_proc(res)
        
    # combine
    combine = Combine()
    proc.add_lc_proc(combine)
    
    # mask
    shower_mask = Mask("Shower")
    proc.add_lc_proc(shower_mask)

    track_mask = Mask("Track")
    proc.add_lc_proc(track_mask)

    # image
    shower_image = Image("Shower")
    proc.add_lc_proc(shower_image)

    track_image = Image("Track")
    proc.add_lc_proc(track_image)

    
