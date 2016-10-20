# Include your header file location
CXXFLAGS += -g -I$(LARLITECV_APPDIR)/ThruMu
CXXFLAGS += -I. $(shell root-config --cflags) -g

CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larlite-config --includes)/../UserDev
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larcv-config --includes)/../app
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)/../app
CXXFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include

# Include your shared object lib location
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g
LDFLAGS += -lprofiler

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = thrumu analyze

all: $(PROGRAMS)

%: %.cxx
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
