import os,sys

if len(sys.argv) != 9:
    print 
    print "PGRAPH_FILE = str(sys.argv[1])"
    print "HIT_FILE    = str(sys.argv[2])"
    print "MCINFO_FILE = str(sys.argv[3])"
    print "OUTPUT_DIR  = str(sys.argv[4])"
    print "REQPDG      = int(sys.argv[5])"
    print "RECLUSTER   = int(sys.argv[6])"
    print "IS_MC       = int(sys.argv[7])"
    print "CONFIG      = str(sys.argv[8])"
    print 
    sys.exit(1)


PGRAPH_FILE = str(sys.argv[1])
HIT_FILE    = str(sys.argv[2])
MCINFO_FILE = str(sys.argv[3])
OUTPUT_DIR  = str(sys.argv[4])
REQPDG      = int(sys.argv[5])
RECLUSTER   = int(sys.argv[6])
IS_MC       = int(sys.argv[7])
CONFIG      = str(sys.argv[8])

num = int(os.path.basename(PGRAPH_FILE).split(".")[0].split("_")[-1])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

#
# Configure Processor
#
proc = llcv.Processor()

#
# handshake
#
dlhs = llcv.DLHandshake()
proc.add_llcv_ana(dlhs)


import reclusterElec
if RECLUSTER == 1:
    #
    # make shower hits
    #
    dlhs = llcv.ShowerHitMaker()
    proc.add_llcv_ana(dlhs)

    #
    # recluster
    #
    # make the raw cluster
    raw_cluster = reclusterElec.ClusterRaw()
    proc.add_ll_ana(raw_cluster)
    
    # correlate raw to dl
    dlraw = reclusterElec.CorrelateDLRaw()
    proc.add_ll_ana(dlraw)

    # recluster seeds
    recluster = reclusterElec.Merge()
    proc.add_ll_ana(recluster)

else:
    #
    # copy dl to dlraw producer
    #
    dlcopy = reclusterElec.CopyDL2DLRaw()
    proc.add_ll_ana(dlcopy)
    
#
# run shower reco
#
from showerRecoDL import DLShowerReco3D
print "...set REQPDG=%d" % REQPDG
dlshr3d = DLShowerReco3D(REQPDG,IS_MC)
proc.add_ll_ana(dlshr3d)

#
# configure
#
proc.configure(CONFIG)
proc.add_lcv_input_file(PGRAPH_FILE)
proc.add_ll_input_file(HIT_FILE)
if MCINFO_FILE != "INVALID":
    proc.add_ll_input_file(MCINFO_FILE)

proc.set_output_ll_name(os.path.join(OUTPUT_DIR,"shower_reco_out_%d.root" % num))
proc.initialize()

#
# run over larcv entries
#
proc.batch_process_lcv()

proc.finalize()

sys.exit(0)
