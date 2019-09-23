{  
  printf("Loading LISE tools...\n");
  gROOT->ProcessLine(".L LISETools/EnergyLossModule.cpp");
  gROOT->ProcessLine(".L LISETools/nuclear_masses.cpp");
  gROOT->ProcessLine(".L LISETools/RelativisticKinematics.cpp");
  
  printf("Loading ROOT Event class...\n");
  gROOT->ProcessLine(".L ./lib/libUNISRootEvent.so");
}
