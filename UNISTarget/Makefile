CFLAGS    :=`root-config --cflags --libs` `nuclust-config --cflags --libs`

DllSuf    := so
SrcSuf    := cpp
ObjSuf    := o

INCLUDES  := -I$(UNISimSrc)
INCLUDES  += -I$(UNISimSrc)UNISShared
INCLUDES  += -I$(UNISimSrc)generator
INCLUDES  += -I$(UNISimSrc)UNISTarget

RLIBS     := -L$(UNISimSrc)lib

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

UNISTARGET_LIB := libUNISTarget.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(UNISTARGET_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(UNISTARGET_LIB): $(OBJS)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS) $(RLIBS)"

install:
	@echo "cp $(UNISTARGET_LIB) $(UNISimSrc)lib/."
	@cp $(UNISTARGET_LIB) $(UNISimSrc)lib/.

.PHONY: distclean
distclean:
	@rm -f $(UNISTARGET_LIB) $(OBJS)
	@rm -f $(UNISimSrc)lib/$(UNISTARGET_LIB)

.PHONY: clean
clean:
	@$(RM) -f $(UNISTARGET_LIB) $(OBJS)
	@$(RM) -f *.expand

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
