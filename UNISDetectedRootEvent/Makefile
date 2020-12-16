CFLAGS    :=`root-config --cflags --libs`

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o
PcmSuf    := _rdict.pcm

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISDetectedRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISDetectedRootEvent_LIB := libUNISDetectedRootEvent.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

UNISDetectedRootEvent_DICT  := UNISDetectedRootEventDict.$(SrcSuf)
UNISDetectedRootEvent_DICTH := $(UNISDetectedRootEvent_DICT:.$(SrcSuf)=.h)
UNISDetectedRootEvent_DICTO := $(UNISDetectedRootEvent_DICT:.$(SrcSuf)=.$(ObjSuf))
UNISDetectedRootEvent_PCM   := UNISDetectedRootEventDict$(PcmSuf)

UNISDetectedRootEvent_HDRS  := $(wildcard UNIS*.$(SrcSuf))
UNISDetectedRootEvent_HDRS  := $(UNISDetectedRootEvent_HDRS:.$(SrcSuf)=.h) UNISDetectedRootEventLinkDef.h
UNISDetectedRootEvent_HDRS  := $(filter-out $(UNISDetectedRootEvent_DICTH),$(UNISDetectedRootEvent_HDRS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(UNISDetectedRootEvent_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf) .$(PcmSuf)

$(UNISDetectedRootEvent_DICT): $(UNISDetectedRootEvent_HDRS)
	@echo "Generating dictionary $@..."
	rootcling -f $@ -p $(INCLUDES) $^

$(UNISDetectedRootEvent_LIB): $(OBJS) $(UNISDetectedRootEvent_DICTO)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"

install:
	@echo "cp $(UNISDetectedRootEvent_LIB) $(UNISDetectedRootEvent_PCM) $(UNISimSrc)lib/."
	@cp $(UNISDetectedRootEvent_LIB) $(UNISDetectedRootEvent_PCM) $(UNISimSrc)lib/.

.PHONY: distclean
distclean:
	@rm -f $(UNISDetectedRootEvent_LIB) $(UNISDetectedRootEvent_DICT) $(UNISDetectedRootEvent_DICTH) $(UNISDetectedRootEvent_PCM) $(OBJS)
	@rm -f $(UNISimSrc)lib/$(UNISDetectedRootEvent_LIB) $(UNISimSrc)lib/$(UNISDetectedRootEvent_PCM) $(OBJS)

.PHONY: clean
clean:
	@$(RM) -f $(UNISDetectedRootEvent_LIB) $(UNISDetectedRootEvent_DICT) $(UNISDetectedRootEvent_DICTH) $(UNISDetectedRootEvent_PCM) $(OBJS)
	@$(RM) -f *.expand

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
