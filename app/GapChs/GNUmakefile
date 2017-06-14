#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = GapChs

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

#HEADERS += IOManager.inl

# include options for this package
INCFLAGS  = -g -I.                       #Include itself
INCFLAGS += -I$(LARCV_BASEDIR)/../
INCFLAGS += -I$(LARLITE_BASEDIR)/../
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += $(shell larlitecv-config --includes)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += $(shell larcv-config --libs)
# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

