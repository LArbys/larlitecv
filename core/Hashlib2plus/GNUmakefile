#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = Hashlib2plus

CXX      =  g++

CFLAGS   = -g -fPIC
LDFLAGS  = -shared
INCFLAGS = -I.

LIB     = ${LARLITECV_LIBDIR}/libLArliteCV${PACKAGE_NAME}.so
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out LinkDef.h, $(wildcard *.h))
OBJS    = $(addprefix ${LARLITECV_BUILDDIR}/${PACKAGE_NAME}/, $(SRCS:.cxx=.o) )

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

#all: $(LIB)
#
#${LARLITECV_BUILDDIR}/${PACKAGE_NAME}/%.o : %.cxx %.h
#	@mkdir -p ${LARLITECV_BUILDDIR}/${PACKAGE_NAME}
#	$(CXX) $(CFLAGS) $(INCFLAGS) -c $< -o $@
#
#$(LIB): $(OBJS)
#	$(CXX) $(LDFLAGS) -o $@ $^

pkg_build: # nothing extra to do in this package

pkg_clean: # nothing extra to do in this package
	@rm $(LIB) ${LARLITECV_BUILDDIR}/${PACKAGE_NAME}/*.o
