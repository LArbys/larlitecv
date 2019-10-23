import os, sys

if len(sys.argv) != 5:
    print 
    print "INFILE      = str(sys.argv[1])"
    print "HIT_FILE    = str(sys.argv[2])"
    print "MCINFO_FILE = str(sys.argv[3])"
    print "OUTDIR      = str(sys.argv[4])"
    print 
    sys.exit(1)

from larlite import larlite as fmwk

INFILE      = str(sys.argv[1])
HIT_FILE    = str(sys.argv[2])
MCINFO_FILE = str(sys.argv[3])
OUTDIR      = str(sys.argv[4])

num = int(INFILE.split(".")[0].split("_")[-1])
OUTFILE = os.path.join(OUTDIR,"showerqualsingle_%d.root" % num)

my_proc = fmwk.ana_processor()
my_proc.add_input_file(INFILE)

if MCINFO_FILE != "INVALID":
    my_proc.add_input_file(MCINFO_FILE)

my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.set_ana_output_file(OUTFILE)
my_proc.set_output_file("")

sq_module = fmwk.ShowerQuality_DL()
sq_module.SetVertexProducer("dlraw")
sq_module.SetShowerProducer("showerreco")
my_proc.add_process(sq_module)

my_proc.run()

sys.exit(0)
