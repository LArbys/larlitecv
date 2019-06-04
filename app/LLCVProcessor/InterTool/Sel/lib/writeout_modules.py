import ROOT
from larlite import larlite as ll
from larcv import larcv as lcv

def InterWriteOutLL():
    name = "InterWriteOutLL"
    return ll.InterWriteOutLL()

def InterWriteOutLC():
    name = "InterWriteOutLC"
    return lcv.InterWriteOutLC(name)

def attach_writeout(proc):
    ll_writeout = InterWriteOutLL()
    proc.add_ll_ana(ll_writeout)

    lc_writeout = InterWriteOutLC()
    proc.add_lc_proc(lc_writeout)

    return
