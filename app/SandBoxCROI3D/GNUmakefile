
# NAME OF LIBRARY/PACKAGE
PACKAGE_NAME = SandBoxCROI3D

# ADD BINDARY NAMES HERE
PROGRAMS = test_croi3d

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
CXXFLAGS += $(shell larcv-config --includes)/../app/ANN
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
ifeq ($(LARCV_OPENCV),1)
CXXFLAGS += -I$(OPENCV_INCDIR) -DUSE_OPENCV
endif
INCFLAGS = $(CXXFLAGS)

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

%: %.cxx $(LARLITECV_BUILDDIR)/lib/liblarlitecv.so
	make --directory=$(LARLITECV_BASEDIR)
	@echo '<< x compiling x' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM

