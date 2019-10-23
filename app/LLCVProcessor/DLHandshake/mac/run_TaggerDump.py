import os, sys

if len(sys.argv) != 4:
    print 
    print "TAGGER_FILE = str(sys.argv[1])"
    print "NUM         = str(sys.argv[2])"
    print "OUTDIR      = str(sys.argv[3])"
    print 
    sys.exit(1)

from larlitecv import larlitecv
from larlite import larlite as fmwk

TAGGER_FILE = str(sys.argv[1])
NUM         = str(sys.argv[2])
OUTDIR      = str(sys.argv[3])

OUTFILE = os.path.join(OUTDIR,"tagger_dump_%s.root" % NUM)

my_proc = fmwk.ana_processor()

my_proc.add_input_file(TAGGER_FILE)

my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.set_ana_output_file(OUTFILE)
my_proc.set_output_file("")

mc_module = fmwk.TaggerDump()
my_proc.add_process(mc_module)

my_proc.run()

sys.exit(0)
