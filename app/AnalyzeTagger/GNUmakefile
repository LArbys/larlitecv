# Include your header file location
CXXFLAGS += -I. $(shell root-config --cflags) -g

CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += -I$(LARLITE_BASEDIR)/../
CXXFLAGS += -I$(LARLITE_BASEDIR)/UserDev
CXXFLAGS += -I$(LARLITE_BASEDIR)/UserDev/SelectionTool
CXXFLAGS += -I$(LARLITE_BASEDIR)/UserDev/BasicTool
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += -I$(LARCV_BASEDIR)/app/Supera/APILArLite
CXXFLAGS += -I$(LARCV_BASEDIR)/app
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += -I$(LARLITECV_BASEDIR)/app

ifeq ($(LARCV_OPENCV),1)
CXXFLAGS += -I$(OPENCV_INCDIR) -DUSE_OPENCV
endif

# Include your shared object lib location
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_FhiclLite
LDFLAGS += $(shell larlite-config --libs) -lSelectionTool_OpT0FinderAna -lSelectionTool_OpT0FinderApp \
	-lSelectionTool_OpT0PhotonLibrary -lSelectionTool_OpT0FinderAlgorithms -lSelectionTool_OpT0FinderBase 
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -lMinuit -g

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = run_pixel_analysis
EXCLUDE_PROGRAMS = run_pixel_analysis_crb_with_double_covered_pixels run_pixel_analysis_for_Taritree  plot_flash_timing
PROGRAM_SRC = $(PROGRAMS:=.cxx)
EX_PROGRAM_SRC = $(EXCLUDE_PROGRAMS:=.cxx)
SOURCES = $(filter-out $(PROGRAM_SRC) $(EX_PROGRAM_SRC), $(wildcard *.cxx))
OBJS    = $(SOURCES:.cxx=.o)

all:		$(PROGRAMS)

%.o: %.cxx
	@echo '<<compiling' $@'>>'
	$(CXX) $(CXXFLAGS) -c $^ -o $@

#$(PROGRAMS): %.cxx $(OBJS)
%: %.cxx $(OBJS)
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f *.o $(PROGRAMS)

