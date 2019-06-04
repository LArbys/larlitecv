CC=g++
CFLAGS=-g
INCFLAGS =-I.
INCFLAGS+=$(shell root-config --cflags)
INCFLAGS+=-I$(OPENCV_INCDIR)
INCFLAGS+=$(shell geo2d-config --includes)
INCFLAGS+=$(shell larlite-config --includes)
INCFLAGS+=-I$(LARLITE_USERDEVDIR)
INCFLAGS+=-I$(LAROPENCV_BASEDIR)
INCFLAGS+=$(shell larcv-config --includes)
INCFLAGS+=$(shell larlitecv-config --includes)

LDFLAGS =-L.

LDFLAGS+=$(shell root-config --libs)
LDFLAGS+=-L$(OPENCV_LIBDIR)
LDFLAGS+=-lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc
LDFLAGS+=$(shell geo2d-config --libs)
LDFLAGS+=$(shell larlitecv-config --libs)
LDFLAGS+=$(shell larlite-config --libs)
LDFLAGS+=-lBasicTool_FhiclLite
LDFLAGS+=-lLArOpenCV_ImageClusterBase
LDFLAGS+=-lLArOpenCV_ImageClusterAlgoFunction
LDFLAGS+=-lLArOpenCV_ImageClusterAlgoData
LDFLAGS+=-lLArOpenCV_ImageClusterAlgoClass
LDFLAGS+=-lLArOpenCV_ImageClusterAlgoModule
LDFLAGS+=-lSelectionTool_OpT0FinderBase
LDFLAGS+=-lSelectionTool_OpT0FinderAlgorithms
LDFLAGS+=-lSelectionTool_OpT0FinderAna
LDFLAGS+=-lSelectionTool_OpT0FinderApp
LDFLAGS+=-lSelectionTool_OpT0PhotonLibrary

SOURCES=tagcosmic.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=tagcosmic

all: $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJECTS): $(SOURCES)
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm *.o tagcosmic tagcosmic.o
