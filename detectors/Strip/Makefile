CFLAGS    :=`root-config --cflags --libs` -lGeom -lGeomPainter -lGeomBuilder -lEve -lGed -lRGL -lEG -lTreePlayer -lMathMore `nuclust-config --cflags --libs`

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISFramework
INCLUDES  += -I$(UNISimSrc)detectors/DetectionSetup
INCLUDES  += -I$(UNISimSrc)detectors/Strip
INCLUDES  += -I$(UNISimSrc)detectors/Lamp
INCLUDES  += -I$(UNISimSrc)UNISRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISSTRIP_LIB := libUNISStrip.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(UNISSTRIP_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(UNISSTRIP_LIB): $(OBJS)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"

install:
	@echo "cp $(UNISSTRIP_LIB) $(UNISimSrc)lib/."
	@cp $(UNISSTRIP_LIB) $(UNISimSrc)lib/.

.PHONY: distclean
distclean:
	@rm -f $(UNISSTRIP_LIB) $(OBJS)
	@rm -f $(UNISimSrc)lib/$(UNISSTRIP_LIB)

.PHONY: clean
clean:
	@$(RM) -f $(UNISSTRIP_LIB) $(OBJS)
	@$(RM) -f *.expand

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
