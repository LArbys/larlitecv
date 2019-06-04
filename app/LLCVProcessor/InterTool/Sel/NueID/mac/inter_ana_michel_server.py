import os,sys,gc

if len(sys.argv) != 9:
    print "SUPERA_FILE = str(sys.argv[1])"
    print "SSNET_FILE  = str(sys.argv[2])"
    print "TAGGER_FILE = str(sys.argv[3])"
    print "VTX_FILE    = str(sys.argv[4])"
    print "RECO2D_FILE = str(sys.argv[5])"
    print "CONFIG_FILE = str(sys.argv[6])"
    print "NUM         = str(sys.argv[7])"
    print "OUT_DIR     = str(sys.argv[8])"
    sys.exit(1)

SUPERA_FILE = str(sys.argv[1])
SSNET_FILE  = str(sys.argv[2])
TAGGER_FILE = str(sys.argv[3])
VTX_FILE    = str(sys.argv[4])
RECO2D_FILE = str(sys.argv[5])
CONFIG_FILE = str(sys.argv[6])
NUM         = str(sys.argv[7])
OUT_DIR     = str(sys.argv[8])

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

driver.SetOutputFilename(os.path.join(OUT_DIR,"michelid_ana_%s.root" % NUM));

selection = llcv.SelMichelID()
selection.SetDebug(False)

driver.AddSelection(selection);

# process
proc.add_llcv_ana(imod)

proc.configure(CONFIG_FILE)

proc.add_lcv_input_file(SUPERA_FILE)
proc.add_lcv_input_file(SSNET_FILE)
proc.add_lcv_input_file(TAGGER_FILE)
proc.add_lcv_input_file(VTX_FILE)

proc.add_ll_input_file(RECO2D_FILE)

proc.set_output_ll_name(os.path.join(OUT_DIR,"michelid_ll_out_%s.root" % NUM))
proc.dataco().set_outputfile(os.path.join(OUT_DIR, "michelid_lcv_out_%s.root" % NUM),"larcv")

proc.initialize()

proc.batch_process_lcv_reverse()

proc.finalize()

sys.exit(0)

