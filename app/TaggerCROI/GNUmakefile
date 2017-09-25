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
#INCFLAGS += $(shell larcv-config --includes)/../app
#INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/GeoAlgo
#INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/FhiclLite
INCFLAGS += $(shell larlitecv-config --includes) 
INCFLAGS += -I$(LARLITECV_BASEDIR)/app # Add this to obtain access to the 'GeneralFlashMatchAlgo' machinery.
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
INCFLAGS += -I$(LARLITE_BASEDIR)/../

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlitecv-config --libs) # Add this to include the analog of everything for larlitecv that is here for larlite/LArCV.
# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

