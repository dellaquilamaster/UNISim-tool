CFLAGS    :=`root-config --cflags --libs` -lGeom -lGeomPainter -lGeomBuilder -lEve -lGed -lRGL -lEG -lTreePlayer -lMathMore `nuclust-config --cflags --libs`

PROG      :=exec_UNISim-tool.exe

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o

INCLUDES  := -I$(UNISimSrc)
INCLUDES  := -I$(UNISimSrc)main
INCLUDES  += -I$(UNISimSrc)UNISFramework
INCLUDES  += -I$(UNISimSrc)detectors/DetectionSetup
INCLUDES  += -I$(UNISimSrc)detectors/Strip
INCLUDES  += -I$(UNISimSrc)detectors/Lamp
INCLUDES  += -I$(UNISimSrc)detectors/Fazia
INCLUDES  += -I$(UNISimSrc)detectors/Silicon
INCLUDES  += -I$(UNISimSrc)detectors/Oscar
INCLUDES  += -I$(UNISimSrc)generator
INCLUDES  += -I$(UNISimSrc)UNISRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared
INCLUDES  += -I$(UNISimSrc)LISETools

RLIBS     := -L$(UNISimSrc)lib -lUNISFramework -lUNISGenerator -lUNISDetectionSetup -lUNISStrip -lUNISLamp -lUNISFaziaQuartet -lUNISSiliconDetector -lUNISOscar -lUNISRootEvent -lUNISShared

SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

.PHONEY : all build
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all : build

build : $(PROG)

$(PROG) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(RLIBS) $(CFLAGS)

install:
	@cp $(PROG) $(UNISimSrc)/$(PROG)

.PHONY: clean
clean :
	rm -f $(OBJS) $(PROG)

.PHONY: distclean
distclean: clean
	@rm -f $(UNISimSrc)/$(PROG)

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
