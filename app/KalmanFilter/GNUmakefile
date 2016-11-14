#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = KalmanFilter

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h 
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

#HEADERS += IOManager.inl

# include options for this package
INCFLAGS  = -g -I.                       #Include itself
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += $(shell larcv-config --includes)/../app
INCFLAGS += $(shell larlite-config --includes) 
INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/GeoAlgo
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(KALMAN_INCDIR)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo $(KALMAN_LIBDIR)/libkalman.a
# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

