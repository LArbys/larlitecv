# Include your header file location
CXXFLAGS += -I. $(shell root-config --cflags)

CXXFLAGS += -I$(LARLITE_BASEDIR)/../
CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larlite-config --includes)/../UserDev
CXXFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool
CXXFLAGS += $(shell larlite-config --includes)/../UserDev/SelectionTool
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larcv-config --includes)/../app
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)/../app
CXXFLAGS += -I$(LARLITECV_BASEDIR)/app # Add this so that 'GeneralFlashMatchAlgo/GeneralFlashMatchAlgoConfig.h'.                                                                           
#CXXFLAGS += -I$(KALMAN_INCDIR)
CXXFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include

ifeq ($(LARCV_OPENCV),1)
CXXFLAGS += -I$(OPENCV_INCDIR) -DUSE_OPENCV
endif

# Include your shared object lib location
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_FhiclLite
LDFLAGS += $(shell larlite-config --libs) -lSelectionTool_OpT0FinderAna -lSelectionTool_OpT0FinderApp \
	-lSelectionTool_OpT0PhotonLibrary -lSelectionTool_OpT0FinderAlgorithms -lSelectionTool_OpT0FinderBase 
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g
#LDFLAGS += -L$(ANN_LIBDIR) -lANN
#LDFLAGS += -lprofiler

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = run_tagger run_cosmictagger

all: $(PROGRAMS)

%: %.cxx
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
