PACKAGE_NAME = InterTool_Sel_FlashMatch

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I.      
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)

# flash matching
INCFLAGS += -I$(LARLITE_USERDEVDIR)/BasicTool
INCFLAGS += -I$(LARLITE_USERDEVDIR)/SelectionTool

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

# flash matching
LDFLAGS += -lSelectionTool_OpT0FinderBase
LDFLAGS += -lSelectionTool_OpT0FinderAlgorithms
LDFLAGS += -lSelectionTool_OpT0FinderAna
LDFLAGS += -lSelectionTool_OpT0FinderApp
LDFLAGS += -lSelectionTool_OpT0PhotonLibrary


# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
pkg_clean:

