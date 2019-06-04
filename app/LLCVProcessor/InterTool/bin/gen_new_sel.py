#!/usr/bin/python
import os, sys


if len(sys.argv) != 2:
    print
    print "\tInterTool Selection template"
    print
    print "\tMust provide:"
    print "\tNAME = str(sys.argv[1])"
    print 
    print "\t\tthank ya bye"
    print
    sys.exit(1)

print
print "\tInterTool Selection template"
print
NAME = str(sys.argv[1])
SELNAME = "Sel%s" % NAME

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)

TEMPLATE_PATH = os.path.join(BASE_PATH,"..","Sel","template")

CLASS_H   = os.path.join(TEMPLATE_PATH,"class_h.template")
CLASS_CXX = os.path.join(TEMPLATE_PATH,"class_cxx.template")
LINKDEF_H = os.path.join(TEMPLATE_PATH,"linkdef_h.template")
MAKEFILE  = os.path.join(TEMPLATE_PATH,"makefile.template")

MKDIR = "mkdir -p %s"

SELDIR = os.path.join(TEMPLATE_PATH,"..")
OUTDIR = os.path.join(SELDIR,NAME)

#
# Make output directory
#
if os.path.exists(OUTDIR) == True:
    print
    print "This directory @ %s already exists" % OUTDIR
    print
    sys.exit(1)

SS = MKDIR % OUTDIR
os.system(SS)


#
# class header and source file
#
with open(os.path.join(OUTDIR,"%s.h" % SELNAME),"w+") as fout:

    data = ""
    with open(CLASS_H,'r') as fin:        
        data = fin.read()

    data = data.replace("AAA",SELNAME)
    data = data.replace("BBB",SELNAME.upper())
    fout.write(data)

with open(os.path.join(OUTDIR,"%s.cxx" % SELNAME),"w+") as fout:

    data = ""
    with open(CLASS_CXX,'r') as fin:        
        data = fin.read()

    data = data.replace("AAA",SELNAME)
    data = data.replace("BBB",SELNAME.upper())
    fout.write(data)
    
#
# linkdef
#   
with open(os.path.join(OUTDIR,"LinkDef.h"),"w+") as fout:

    data = ""
    with open(LINKDEF_H,'r') as fin:        
        data = fin.read()

    data = data.replace("AAA",SELNAME)
    fout.write(data)

#
# GNUMakefile
#
with open(os.path.join(OUTDIR,"GNUmakefile"),"w+") as fout:

    data = ""
    with open(MAKEFILE,'r') as fin: 
        data = fin.read()

    data = data.replace("AAA",SELNAME)
    fout.write(data)


#
# Add me to the Sel toplevel makefile
#
with open(os.path.join(SELDIR,"GNUmakefile"),'r') as fin:
    data = fin.read()

data = data.replace("Toy","Toy %s" % NAME)

with open(os.path.join(SELDIR,"GNUmakefile"),'w+') as fout:
    fout.write(data)

print "\t...DONE!"
print
print "\t\tthank ya bye"
print
sys.exit(0)
