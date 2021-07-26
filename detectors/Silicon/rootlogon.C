{
  gROOT->ProcessLine(".L ../DetectionSetup/UNISDetectionSetup.cpp"); 
  gROOT->ProcessLine(".L ../DetectionSetup/UNISDetectionUnit.cpp");
  gROOT->ProcessLine(".L ../Oscar/UNISOscarTelescope.cpp");
  gROOT->ProcessLine(".L ../Strip/UNISStripSingleSidedDetector.cpp");
  gROOT->ProcessLine(".L UNISSiliconPhotoDiode.cpp");
}
