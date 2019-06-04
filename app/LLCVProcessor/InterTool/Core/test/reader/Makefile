CC=g++

INCFLAGS = -I.
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell geo2d-config --includes)
INCFLAGS += -I$(LAROPENCV_BASEDIR)
INCFLAGS += -I$(LARCV_BASEDIR)/app/LArOpenCVHandle
INCFLAGS += $(shell root-config --cflags)
INCFLAGS += $(shell larlitecv-config --includes)

LDFLAGS = -L.
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += -lBasicTool_FhiclLite
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell geo2d-config --libs)
LDFLAGS += -lLArOpenCV_Core
LDFLAGS += -lLArOpenCV_ImageClusterBase
LDFLAGS += -lLArOpenCV_ImageClusterAlgoFunction
LDFLAGS += -lLArOpenCV_ImageClusterAlgoData
LDFLAGS += -lLArOpenCV_ImageClusterAlgoClass
LDFLAGS += $(shell root-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += -lTreePlayer

SOURCES=test.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=test

all: $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJECTS): $(SOURCES)
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm *.o test
