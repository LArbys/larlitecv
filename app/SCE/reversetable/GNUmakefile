# Include your header file location
CXXFLAGS += -I$(LARLITECV_APPDIR)/SCE
CXXFLAGS += -I. $(shell root-config --cflags)

CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)

# Include your shared object lib location
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g
LDFLAGS += -L$(ANN_LIBDIR) -lANN

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = makeReverseMapData

all: $(PROGRAMS)

%: %.cxx
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
