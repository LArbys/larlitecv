#dirtyness by vic
CXX=g++
CXXFLAGS += -O3 #-g

######################
#############INCLUDES#
######################

##ROOT
BASIC_ROOTINC=$(shell root-config --cflags)
CXXINCS +=$(BASIC_ROOTINC)

##larlite
LARLITE_INCDIR=$(shell larlite-config --includes)
CXXINCS += $(LARLITE_INCDIR)

CXXINCS += -I$(LARLITE_USERDEVDIR)/BasicTool
CXXINCS += -I$(LARLITE_USERDEVDIR)/SelectionTool

##larcv
LARCV_INCDIR=$(shell larcv-config --includes)
CXXINCS += $(LARCV_INCDIR)

##larlitecv
LARLITECV_INCDIR=$(shell larlitecv-config --includes)
CXXINCS += $(LARLITECV_INCDIR)

##geo2d
CXXINCS += -I$(GEO2D_BASEDIR)

##laropencv
CXXINCS += -I$(LAROPENCV_BASEDIR)

##opencv
CXXINCS +=-I$(OPENCV_INCDIR)

#######################
#############LIBRARIES#
#######################

##ROOT
BASIC_ROOTLIBS=$(shell root-config --libs)
CXXLIBS += $(BASIC_ROOTLIBS)

##larlite
CXXLIBS += $(shell larlite-config --libs)

##larcv
CXXLIBS += $(shell larcv-config --libs)

##larlitecv
CXXLIBS += $(shell larlitecv-config --libs)

CXXLIBS += -lSelectionTool_OpT0FinderBase
CXXLIBS += -lSelectionTool_OpT0FinderAlgorithms
CXXLIBS += -lSelectionTool_OpT0FinderAna
CXXLIBS += -lSelectionTool_OpT0FinderApp
CXXLIBS += -lSelectionTool_OpT0PhotonLibrary

##laropencv
CXXLIBS += -L$(LARLITE_LIBDIR)
CXXLIBS += -lLArOpenCV_Core
CXXLIBS += -lLArOpenCV_ImageClusterBase
CXXLIBS += -lLArOpenCV_ImageClusterAlgoFunction
CXXLIBS += -lLArOpenCV_ImageClusterAlgoData
CXXLIBS += -lLArOpenCV_ImageClusterAlgoClass

##geo2d
CXXLIBS += -L$(GEO2D_LIBDIR)
CXXLIBS += -lGeo2D_Core

##opencv
CXXLIBS +=-L$(OPENCV_LIBDIR)
OPENCV_LIBS=-lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc
CXXLIBS += $(OPENCV_LIBS)

#############OUTPUT
SWFILE=test.cxx
OWFILE=test.o
BWFILE=test

all:
	$(CXX) $(CXXFLAGS) $(CXXINCS) -c $(SWFILE) -o $(OWFILE)
	$(CXX) -o $(BWFILE) $(OWFILE) $(CXXLIBS)
clean:
	rm -rf $(OWFILE) $(BWFILE)
	rm -rf *.png
