#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = DLHandshake

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I.      
INCFLAGS += -I$(GEO2D_BASEDIR)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LAROPENCV_BASEDIR)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)

LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += -lLArOpenCV_Core
LDFLAGS += -lLArOpenCV_ImageClusterBase
LDFLAGS += $(shell larcv-config --libs)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

