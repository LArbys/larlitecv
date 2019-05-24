#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = StopMu

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

#HEADERS += IOManager.inl

# include options for this package
INCFLAGS  = -g -I./
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += $(shell larcv-config --includes)/../app
INCFLAGS += $(shell larlite-config --includes) 
#INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/GeoAlgo
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include

ifeq ($(LARCV_OPENCV),1)
INCFLAGS += -I$(OPENCV_INCDIR) -DUSE_OPENCV
endif



# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

ifeq ($(LARCV_OPENCV),1)
LDFLAGS += -L$(OPENCV_LIBDIR) -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc
endif
LDFLAGS += $(shell larcv-config --libs)
#LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

