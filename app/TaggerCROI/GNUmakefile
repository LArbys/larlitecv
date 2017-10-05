#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = TaggerCROI

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

#HEADERS += IOManager.inl

# include options for this package
INCFLAGS  = -g -I.                       #Include itself
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell basictool-config --includes)
INCFLAGS += $(shell seltool-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
INCFLAGS += -I$(LARLITE_BASEDIR)/../
INCFLAGS += -I$(LARLITECV_BASEDIR)/app
# needs opencv
INCFLAGS += -DUSE_OPENCV=1
INCFLAGS += -I$(LAROPENCV_BASEDIR)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

