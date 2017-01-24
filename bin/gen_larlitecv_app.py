#!/usr/bin/env python
import os,sys

# makes minimial app
if len(sys.argv)!=2:
    print "usage: make_new_larlitecv_app.py [app name]"
    sys.exit(-1)
    
appname = sys.argv[1]

print "Will generate an app in the app folder with name: ",appname
print "[Enter to continue]"
raw_input()


minimial_gnumake="""
# NAME OF LIBRARY/PACKAGE
PACKAGE_NAME = %s

# ADD BINDARY NAMES HERE
PROGRAMS = run_%s

PROGRAM_SOURCE = $(addsuffix .cxx, $(PROGRAMS))
SOURCES = $(filter-out $(PROGRAM_SOURCE), $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# Include your header file location
# and include compiler options for this package
CXXFLAGS += -I. $(shell root-config --cflags) -g
CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larlite-config --includes)/../UserDev
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larcv-config --includes)/../app
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)/../app

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Include your shared object lib location
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.

all: $(PROGRAMS)

$(PROGRAMS):
	make --directory=../..
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
"""

minimal_linkdef="""
//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace larlitecv;
#pragma link C++ namespace larlitecv::%s;
//ADD_NEW_CLASS ... do not change this line
#endif
"""

minimal_class_header="""
#ifndef __larlitecv_%s_h__
#define __larlitecv_%s_h__

namespace larlitecv {

  class %s {
  public:
    %s() {};
    virtual ~%s() {};

  };

}
#endif
"""

minimal_class_source="""
#include "%s.h"

namespace larlitecv {

}
"""

minimal_main="""
#include <iostream>

#include "%s.h"

int main( int nargs, char** argv ) {

  return 0;
}
"""

appdir = os.environ["LARLITECV_BASEDIR"]+"/app/%s"%(appname)
os.system("mkdir -p %s"%(appdir))

header = file(appdir+"/%s.h"%(appname),'w')
source = file(appdir+"/%s.cxx"%(appname),'w')
main   = file(appdir+"/run_%s.cxx"%(appname),'w')

header_tuple = (appname,)*5
print header_tuple

print >> header,minimal_class_header%header_tuple
print >> source,minimal_class_source%(appname)
print >> main,minimal_main%(appname)


gnumakefile = minimial_gnumake % ( appname, appname )

linkdef = file(appdir+"/LinkDef.h",'w')
print >>linkdef,minimal_linkdef%(appname)

makefile = file(appdir+"/GNUmakefile",'w')
print >>makefile,gnumakefile
