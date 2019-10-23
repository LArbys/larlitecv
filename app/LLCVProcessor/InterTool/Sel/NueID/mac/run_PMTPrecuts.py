import os, sys

if len(sys.argv) != 5:
    print 
    print "OPRECO_FILE = str(sys.argv[1])"
    print "CFG         = str(sys.argv[2])"
    print "NUM         = str(sys.argv[3])"
    print "OUTDIR      = str(sys.argv[4])"
    print 
    sys.exit(1)

from larlitecv import larlitecv
from larlite import larlite as fmwk
import ROOT
from ROOT import fcllite
fcllite.PSet

OPRECO_FILE = str(sys.argv[1])
CFG         = str(sys.argv[2])
NUM         = str(sys.argv[3])
OUTDIR      = str(sys.argv[4])

OUTFILE = os.path.join(OUTDIR,"pmt_precut_dump_%s.root" % NUM)

my_proc = fmwk.ana_processor()

my_proc.add_input_file(OPRECO_FILE)

my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.set_ana_output_file(OUTFILE)
my_proc.set_output_file("")

mc_module = fmwk.LEEPreCut()
mc_module.configure(fcllite.CreatePSetFromFile(ROOT.std.string(CFG)))
my_proc.add_process(mc_module)

my_proc.run()

sys.exit(0)
