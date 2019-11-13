CFLAGS    :=`root-config --cflags --libs` -lGeom -lEve -lGed -lRGL -lEG -lTreePlayer -lMathMore

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISFramework
INCLUDES  += -I$(UNISimSrc)detectors/DetectionSetup
INCLUDES  += -I$(UNISimSrc)detectors/Strip
INCLUDES  += -I$(UNISimSrc)detectors/Lamp
INCLUDES  += -I$(UNISimSrc)generator
INCLUDES  += -I$(UNISimSrc)UNISRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared
INCLUDES  += -I$(UNISimSrc)LISETools

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISFRAMEWORK_LIB := libUNISFramework.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(UNISFRAMEWORK_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(UNISFRAMEWORK_LIB): $(OBJS)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"

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