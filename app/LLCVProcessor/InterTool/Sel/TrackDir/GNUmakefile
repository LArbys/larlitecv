PACKAGE_NAME = InterTool_SelTrackDir

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I.      
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)

INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell geo2d-config --includes)
INCFLAGS += -I$(LAROPENCV_BASEDIR)
INCFLAGS += -I$(LARCV_BASEDIR)/app/LArOpenCVHandle
INCFLAGS += $(shell larlitecv-config --includes)

LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += -lBasicTool_FhiclLite
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell geo2d-config --libs)
LDFLAGS += -lLArOpenCV_Core
LDFLAGS += -lLArOpenCV_ImageClusterBase
LDFLAGS += -lLArOpenCV_ImageClusterAlgoFunction
LDFLAGS += -lLArOpenCV_ImageClusterAlgoData
LDFLAGS += -lLArOpenCV_ImageClusterAlgoClass

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

