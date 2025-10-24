CFLAGS    :=`root-config --cflags --libs` -lGeom -lGeomPainter -lGeomBuilder -lEve -lGed -lRGL -lEG -lTreePlayer -lMathMore `nuclust-config --cflags --libs`

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISFramework
INCLUDES  += -I$(UNISimSrc)detectors/DetectionSetup
INCLUDES  += -I$(UNISimSrc)detectors/Strip
INCLUDES  += -I$(UNISimSrc)detectors/Fazia
INCLUDES  += -I$(UNISimSrc)detectors/Lamp
INCLUDES  += -I$(UNISimSrc)detectors/Silicon
INCLUDES  += -I$(UNISimSrc)detectors/Oscar
INCLUDES  += -I$(UNISimSrc)generator
INCLUDES  += -I$(UNISimSrc)UNISRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared
INCLUDES  += -I$(UNISimSrc)UNISTarget

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISFRAMEWORK_LIB := libUNISFramework.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

.$(SrcSuf).$(ObjSuf):
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)

all: $(UNISFRAMEWORK_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(UNISFRAMEWORK_LIB): $(OBJS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)

install:
	@echo "cp $(UNISFRAMEWORK_LIB) $(UNISimSrc)lib/."
	@cp $(UNISFRAMEWORK_LIB) $(UNISimSrc)lib/.

.PHONY: distclean
distclean:
	@rm -f $(UNISFRAMEWORK_LIB) $(OBJS)
	@rm -f $(UNISimSrc)lib/$(UNISFRAMEWORK_LIB)

.PHONY: clean
clean:
	@$(RM) -f $(UNISFRAMEWORK_LIB) $(OBJS)
	@$(RM) -f *.expand

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
