CXX      =  g++

CFLAGS   = -g -fPIC
LDFLAGS  = -shared
INCLUDES = -I.

MODULE  = Hashlib2plus
LIB     = ${LARLITECV_LIBDIR}/libLArliteCV${MODULE}.so
SRCS    = $(wildcard *.cxx)
OBJS    = $(addprefix ${LARLITECV_BUILDDIR}/${MODULE}/, $(SRCS:.cxx=.o) )
HEADERS = $(wildcard *.h)

all: $(LIB)

${LARLITECV_BUILDDIR}/${MODULE}/%.o : %.cxx %.h
	@mkdir -p ${LARLITECV_BUILDDIR}/${MODULE}
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(LIB): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	@rm $(LIB) ${LARLITECV_BUILDDIR}/${MODULE}/*.o
