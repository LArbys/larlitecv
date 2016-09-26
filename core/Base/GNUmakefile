#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = Base

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
#IMPS = $(wildcard *.inl)

#USE_OPENCV=$(LARLITECV_OPENCV)

# include options for this package
INCFLAGS  = -I. $(shell larlitecv-config --includes)  #Include itself

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build: # nothing extra to do in this package
pkg_clean: # nothing extra to do in this package

