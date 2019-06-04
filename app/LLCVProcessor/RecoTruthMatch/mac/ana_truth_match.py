import os,sys

if len (sys.argv) != 6:
    print 
    print 
    print "CONFIG_FILE     = str(sys.argv[1])"
    print "TYPE            = str(sys.argv[2])"
    print "SSNET_FILE      = str(sys.argv[3])"
    print "[TRK|SHR]_FILE  = str(sys.argv[4])"
    print "OUTPUT_DIR      = str(sys.argv[5])"
    print 
    sys.exit(1)

from match_factory import MatchFactory


CONFIG_FILE   = str(sys.argv[1])
TYPE          = str(sys.argv[2])
SSNET_FILE    = str(sys.argv[3])
TRK_SHR_FILE  = str(sys.argv[4])
OUTPUT_DIR    = str(sys.argv[5])

num = int(os.path.basename(TRK_SHR_FILE).split(".")[0].split("_")[-1])

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

mf = MatchFactory(TYPE)
match_algo = mf.get()

proc.add_llcv_ana(match_algo)


proc.configure(CONFIG_FILE)
proc.dataco().set_outputfile(os.path.join(OUTPUT_DIR, mf.output_file(num)),"larcv")
proc.add_lcv_input_file(SSNET_FILE)
proc.add_ll_input_file(TRK_SHR_FILE)
proc.initialize()

#
# run over larcv entries
#
proc.batch_process_lcv()

proc.finalize()

sys.exit(0)
