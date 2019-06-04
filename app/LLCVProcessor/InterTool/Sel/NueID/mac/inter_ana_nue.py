import os,sys,gc

if len(sys.argv) != 8:
    print "SSNET_FILE  = str(sys.argv[1])"
    print "VTX_FILE    = str(sys.argv[2])"
    print "RECO2D_FILE = str(sys.argv[3])"
    print "CONFIG_FILE = str(sys.argv[4])"
    print "NUM         = str(sys.argv[5])"
    print "IS_MC       = int(sys.argv[6])"
    print "OUT_DIR     = str(sys.argv[7])"
    sys.exit(1)

SSNET_FILE  = str(sys.argv[1])
VTX_FILE    = str(sys.argv[2])
RECO2D_FILE = str(sys.argv[3])
CONFIG_FILE = str(sys.argv[4])
NUM         = str(sys.argv[5])
IS_MC       = int(str(sys.argv[6]))
OUT_DIR     = str(sys.argv[7])

import ROOT
from larlitecv import larlitecv
from ROOT import llcv

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
print "Base path at: ",BASE_PATH
sys.path.insert(0,BASE_PATH)

TOP_DIR = os.environ['LARLITECV_BASEDIR']
TOP_DIR = os.path.join(TOP_DIR,"app","LLCVProcessor","InterTool")
sys.path.insert(0,os.path.join(TOP_DIR,"Sel"))

from lib.ssnet_modules import attach_ssnet
from lib.dead_modules import attach_dead
from lib.writeout_modules import attach_writeout

# proc
proc = llcv.Processor()

# attach SSNet
attach_ssnet(proc)

# attach dead wire maker
attach_dead(proc)

# attach the larlite writeout
attach_writeout(proc)

# intermodule
imod = llcv.InterModule()

# configure the driver
driver = imod.Driver()

driver.SetOutputFilename(os.path.join(OUT_DIR,"nueid_ana_%s.root" % NUM));

selection = llcv.SelNueID()
selection.SetIsMC(IS_MC)
selection.SetDebug(True)

driver.AddSelection(selection);

# process
proc.add_llcv_ana(imod)

proc.configure(CONFIG_FILE)

proc.add_lcv_input_file(SSNET_FILE)
proc.add_lcv_input_file(VTX_FILE)
proc.add_ll_input_file(RECO2D_FILE)

proc.dataco().set_outputfile(os.path.join(OUT_DIR, "nueid_lcv_out_%s.root" % NUM),"larcv")
proc.set_output_ll_name(os.path.join(OUT_DIR,"nueid_ll_out_%s.root" % NUM))

proc.initialize()

#proc.batch_process_lcv_reverse(16,1)
#proc.batch_process_lcv_reverse(0,1)
proc.batch_process_lcv_reverse()

proc.finalize()

sys.exit(0)

