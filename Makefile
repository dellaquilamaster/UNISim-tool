CFLAGS    :=`root-config --cflags --libs` -lGeom -lEve -lGed -lRGL -lEG -lTreePlayer -lMathMore

SrcSuf    := cpp
ObjSuf    := o
ExeSuf    := exe
DepSuf    := h
DllSuf    := so
PcmSuf    := _rdict.pcm

ObjDir    := obj/
LibDir    := lib/
SrcDir    := src/
DepDir    := include/

PROG      := $(wildcard exec_*.$(SrcSuf))
PROG      := $(patsubst %.$(SrcSuf), %, $(PROG))

UNISROOTEVENT_DICT  := UNISRootEventDict.$(SrcSuf)
UNISROOTEVENT_DICTH := UNISRootEventDict.$(DepSuf)
UNISROOTEVENT_DICTO := UNISRootEventDict.$(ObjSuf)
UNISROOTEVENT_PCM   := UNISRootEventDict$(PcmSuf)
UNISROOTEVENT_HDRS  := UNISRootEventLinkDef.$(DepSuf)

OBJS      := $(ObjDir)TDetectionSetup.$(ObjSuf) $(ObjDir)TDetectionUnit.$(ObjSuf) $(ObjDir)TStripDetector.$(ObjSuf) $(ObjDir)TLampWedgeDetector.$(ObjSuf) $(ObjDir)EnergyLossModule.$(ObjSuf) $(ObjDir)nuclear_masses.$(ObjSuf) $(ObjDir)RelativisticKinematics.$(ObjSuf) $(ObjDir)UNISRootEvent.$(ObjSuf) $(ObjDir)UNISEventGenerator.$(ObjSuf) $(ObjDir)UNISSequentialDecay.$(ObjSuf) $(ObjDir)UNISRutherfordScattering.$(ObjSuf) $(ObjDir)UNISFramework.$(ObjSuf) $(ObjDir)$(UNISROOTEVENT_DICTO) $(ObjDir)UNISLogo.$(ObjSuf) $(ObjDir)shared.$(ObjSuf)

UNIS_LIB  := libUNISim.$(DllSuf)

RLIBS     := -L./$(LibDir) -lUNISim

DEPS      := $(_OBJS:.$(ObjSuf)=.$(DepSuf))

INCLUDES  := -I./detectors/DetectionSetup/
INCLUDES  += -I./detectors/Strip/
INCLUDES  += -I./detectors/Lamp/
INCLUDES  += -I./LISETools/
INCLUDES  += -I./include/
INCLUDES  += -I./generator/

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3

all: $(PROG) $(UNIS_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(PROG): $(OBJS)
		$(CXX) $(CXXFLAGS) -o ${@}.$(ExeSuf) ${@}.cpp $^ $(SYSLIB) $(CFLAGS)
		
$(UNIS_LIB): $(OBJS)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS)"	
		
$(ObjDir)TDetectionSetup.$(ObjSuf): ./detectors/DetectionSetup/TDetectionSetup.cpp  ./detectors/DetectionSetup/TDetectionSetup.h
	$(CXX) $(CXXFLAGS) -c  ./detectors/DetectionSetup/TDetectionSetup.cpp -o $(ObjDir)TDetectionSetup.$(ObjSuf) $(SYSLIB) $(CFLAGS) 	
	
$(ObjDir)TDetectionUnit.$(ObjSuf): ./detectors/DetectionSetup/TDetectionUnit.cpp  ./detectors/DetectionSetup/TDetectionUnit.h
	$(CXX) $(CXXFLAGS) -c  ./detectors/DetectionSetup/TDetectionUnit.cpp -o $(ObjDir)TDetectionUnit.$(ObjSuf) $(SYSLIB) $(CFLAGS) 	
	
$(ObjDir)TStripDetector.$(ObjSuf): ./detectors/Strip/TStripDetector.cpp  ./detectors/Strip/TStripDetector.h
	$(CXX) $(CXXFLAGS) -c  ./detectors/Strip/TStripDetector.cpp -o $(ObjDir)TStripDetector.$(ObjSuf) $(SYSLIB) $(CFLAGS) 	
	
$(ObjDir)TLampWedgeDetector.$(ObjSuf): ./detectors/Lamp/TLampWedgeDetector.cpp  ./detectors/Lamp/TLampWedgeDetector.h
	$(CXX) $(CXXFLAGS) -c  ./detectors/Lamp/TLampWedgeDetector.cpp -o $(ObjDir)TLampWedgeDetector.$(ObjSuf) $(SYSLIB) $(CFLAGS) 	

$(ObjDir)EnergyLossModule.$(ObjSuf):  ./LISETools/EnergyLossModule.cpp  ./LISETools/EnergyLossModule.h
	$(CXX) $(CXXFLAGS) -c ./LISETools/EnergyLossModule.cpp -o $(ObjDir)EnergyLossModule.$(ObjSuf) $(SYSLIB) $(CFLAGS) 

$(ObjDir)nuclear_masses.$(ObjSuf):  ./LISETools/nuclear_masses.cpp  ./LISETools/nuclear_masses.h
	$(CXX) $(CXXFLAGS) -c ./LISETools/nuclear_masses.cpp -o $(ObjDir)nuclear_masses.$(ObjSuf) $(SYSLIB) $(CFLAGS) 

$(ObjDir)RelativisticKinematics.$(ObjSuf):  ./LISETools/RelativisticKinematics.cpp  ./LISETools/RelativisticKinematics.h
	$(CXX) $(CXXFLAGS) -c ./LISETools/RelativisticKinematics.cpp -o $(ObjDir)RelativisticKinematics.$(ObjSuf) $(SYSLIB) $(CFLAGS) 

$(ObjDir)UNISFramework.$(ObjSuf):  ./$(SrcDir)UNISFramework.cpp  ./$(DepDir)UNISFramework.h
	$(CXX) $(CXXFLAGS) -c ./$(SrcDir)UNISFramework.cpp -o $(ObjDir)UNISFramework.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)shared.$(ObjSuf):  ./$(SrcDir)shared.cpp  ./$(DepDir)shared.h
	$(CXX) $(CXXFLAGS) -c ./$(SrcDir)shared.cpp -o $(ObjDir)shared.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)UNISLogo.$(ObjSuf):  ./$(SrcDir)UNISLogo.cpp  ./$(DepDir)UNISLogo.h
	$(CXX) $(CXXFLAGS) -c ./$(SrcDir)UNISLogo.cpp -o $(ObjDir)UNISLogo.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)UNISEventGenerator.$(ObjSuf):  ./generator/UNISEventGenerator.cpp  ./generator/UNISEventGenerator.h
	$(CXX) $(CXXFLAGS) -c ./generator/UNISEventGenerator.cpp -o $(ObjDir)UNISEventGenerator.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)UNISSequentialDecay.$(ObjSuf):  ./generator/UNISSequentialDecay.cpp  ./generator/UNISSequentialDecay.h
	$(CXX) $(CXXFLAGS) -c ./generator/UNISSequentialDecay.cpp -o $(ObjDir)UNISSequentialDecay.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)UNISRutherfordScattering.$(ObjSuf):  ./generator/UNISRutherfordScattering.cpp  ./generator/UNISRutherfordScattering.h
	$(CXX) $(CXXFLAGS) -c ./generator/UNISRutherfordScattering.cpp -o $(ObjDir)UNISRutherfordScattering.$(ObjSuf) $(SYSLIB) $(CFLAGS) 

$(ObjDir)UNISRootEvent.$(ObjSuf):  ./$(SrcDir)UNISRootEvent.cpp  ./$(DepDir)UNISRootEvent.h ./$(ObjDir)$(UNISROOTEVENT_DICTO)
	$(CXX) $(CXXFLAGS) -c ./$(SrcDir)UNISRootEvent.cpp -o $(ObjDir)UNISRootEvent.$(ObjSuf) $(SYSLIB) $(CFLAGS) 
	
$(ObjDir)$(UNISROOTEVENT_DICTO):  ./$(UNISROOTEVENT_DICT)
	$(CXX) $(CXXFLAGS) -c ./$(UNISROOTEVENT_DICT) -o ./$(ObjDir)$(UNISROOTEVENT_DICTO) $(SYSLIB) $(CFLAGS) 
	
$(UNISROOTEVENT_DICT): ./$(DepDir)UNISRootEvent.h ./$(DepDir)$(UNISROOTEVENT_HDRS)
	rootcling -f $@ -p $(INCLUDES) $^
	
install:
	@cp $(UNIS_LIB) $(LibDir)$(UNIS_LIB)	
	@cp $(UNISROOTEVENT_PCM) $(LibDir)$(UNISROOTEVENT_PCM)	
	
.PHONY: clean
clean:
	@$(RM) *.$(ExeSuf)
	@$(RM) $(UNIS_LIB)
	@$(RM) $(ObjDir)*.$(ObjSuf)
	@$(RM) $(UNISROOTEVENT_DICT) $(UNISROOTEVENT_DICTH) $(UNISROOTEVENT_DICTO) $(UNISROOTEVENT_PCM)

.PHONY: distclean
distclean: clean
	@rm -f $(LibDir)/$(UNIS_LIB)
	@rm -f $(LibDir)/$(UNISROOTEVENT_PCM)
