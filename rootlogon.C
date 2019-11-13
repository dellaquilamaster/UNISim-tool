{  
  printf("Loading LISE tools...\n");
  gROOT->ProcessLine(".L LISETools/EnergyLossModule.cpp");
  gROOT->ProcessLine(".L LISETools/nuclear_masses.cpp");
  gROOT->ProcessLine(".L LISETools/RelativisticKinematics.cpp");
  
  printf("Loading Shared libraries...\n");
  gROOT->ProcessLine(".L ./lib/libLISETools.so");
  gROOT->ProcessLine(".L ./lib/libUNISDetectionSetup.so");
  gROOT->ProcessLine(".L ./lib/libUNISFramework.so");
  gROOT->ProcessLine(".L ./lib/libUNISGenerator.so");
  gROOT->ProcessLine(".L ./lib/libUNISLamp.so");
  gROOT->ProcessLine(".L ./lib/libUNISRootEvent.so");
  gROOT->ProcessLine(".L ./lib/libUNISShared.so");
  gROOT->ProcessLine(".L ./lib/libUNISStrip.so");  
}
