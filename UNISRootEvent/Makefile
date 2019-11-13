CFLAGS    :=`root-config --cflags --libs`

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o
PcmSuf    := _rdict.pcm

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISRootEvent
INCLUDES  += -I$(UNISimSrc)UNISShared

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISRootEvent_LIB := libUNISRootEvent.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

UNISRootEvent_DICT  := UNISRootEventDict.$(SrcSuf)
UNISRootEvent_DICTH := $(UNISRootEvent_DICT:.$(SrcSuf)=.h)
UNISRootEvent_DICTO := $(UNISRootEvent_DICT:.$(SrcSuf)=.$(ObjSuf))
UNISRootEvent_PCM   := UNISRootEventDict$(PcmSuf)

UNISRootEvent_HDRS  := $(wildcard UNIS*.$(SrcSuf))
UNISRootEvent_HDRS  := $(UNISRootEvent_HDRS:.$(SrcSuf)=.h) UNISRootEventLinkDef.h
UNISRootEvent_HDRS  := $(filter-out $(UNISRootEvent_DICTH),$(UNISRootEvent_HDRS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(UNISRootEvent_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf) .$(PcmSuf)

$(UNISRootEvent_DICT): $(UNISRootEvent_HDRS)
	@echo "Generating dictionary $@..."
	rootcling -f $@ -p $(INCLUDES) $^

$(UNISRootEvent_LIB): $(OBJS) $(UNISRootEvent_DICTO)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"

install:
	@echo "cp $(UNISRootEvent_LIB) $(UNISRootEvent_PCM) $(UNISimSrc)lib/."
	@cp $(UNISRootEvent_LIB) $(UNISRootEvent_PCM) $(UNISimSrc)lib/.

.PHONY: distclean
distclean:
	@rm -f $(UNISRootEvent_LIB) $(UNISRootEvent_DICT) $(UNISRootEvent_DICTH) $(UNISRootEvent_PCM) $(OBJS)
	@rm -f $(UNISimSrc)lib/$(UNISRootEvent_LIB) $(UNISimSrc)lib/$(UNISRootEvent_PCM) $(OBJS)

.PHONY: clean
clean:
	@$(RM) -f $(UNISRootEvent_LIB) $(UNISRootEvent_DICT) $(UNISRootEvent_DICTH) $(UNISRootEvent_PCM) $(OBJS)
	@$(RM) -f *.expand

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
