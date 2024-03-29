/* THE UNISim-tool Framework
 * Created by: Daniele Dell'Aquila
 * e-mail:     daniele.dellaquila@irb.hr
 * 
 * The framework is a Unified Simulation tool for simulation of 
 * nuclear physics reactions with Monte Carlo tools.
 * 
 * 
 * 10/06/2019
 * v1.0.beta
 * 03/10/2019
 * v1.1
 * new features:
 *  - energy loss in the target
 *  - possibility to specify output file name
 * 
 * 
 */

#ifndef UNISFRAMEWORK_H
#define UNISFRAMEWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctime>
#include <math.h>
#include <vector>
#include <map>
#include <TRint.h>
#include <TGeoManager.h>
#include <TGeoTrack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <UNISEventGenerator.h>
#include <UNISSequentialDecay.h>
#include <UNISSequentialDecayTwoBody.h>
#include <UNISRutherfordScattering.h>
#include <UNISRootEvent.h>
#include <nuclear_masses.h>
#include <EnergyLossModule.h>
#include <RelativisticKinematics.h>
#include <UNISDetectionSetup.h>
#include <UNISDetectionUnit.h>
#include <UNISStripCluster.h>
#include <UNISStripDetector.h>
#include <UNISFaziaQuartet.h>
#include <UNISFaziaBlock.h>
#include <UNISLampWedgeDetector.h>
#include <UNISLampWedgeMMMDetector.h>
#include <UNISS3Detector.h>
#include <UNISSiliconPhotoDiode.h>
#include <UNISStripSingleSidedDetector.h>
#include <UNISOscarTelescope.h>
#include <UNISCollimatedSilicon.h>
#include <UNISIon.h>
#include <UNISShared.h>

class UNISFramework
{
public :
  UNISFramework(); //Constructor
  UNISFramework(const char *); //Constructor
  ~UNISFramework(); //Destructor
  
  //
  // Configuration methods
  int ConfigureFramework(); //Lauch framework configuration
  void InitTree(); //Initialize output tree
  void SetIterations(int); //Set the number of iterations
  void SetVerbose(bool opt=true); //Set verbose mode
  void SetGraphics(bool opt=true); //Enables openGL graphics
  //
  
  //
  //Methods to run the framework
  int ReadInput(int, char **); //Read input parameters from string
  void ProcessIterations(); //Process fNumEvents iterations of the framework
  void RegisterEvent(std::vector<UNISIon> &); //Fill the current event data structure
  void RegisterEvent(); //Fill the current event data structure from input file
  void FillEvent(); //Fill the tree for the current event
  void EndProcess(); //End iteration process, save tree to file and close file
  //
  
  //
  //Chech the configuration
  void PrintConfiguration() const;
  //
  
private :
  //
  //Framework Data Members
  const std::string fConfigurationFile;
  std::string fOutputFolder;
  std::string fOutputFileName;
  std::string fInputFileName;
  bool fGenerateData;
  bool fGraphics;
  bool fAdvancedGraphics;
  bool fVerbose;
  double fStartTime;
  //
  
  //
  //Beam Data
  UNISIon fTheBeam;
  double fBeamEnergy; //MeV Total kinetic energy
  double fBeamEnergyMidTarget; //MeV Total kinetic energy at mid target (useful for the reaction)
  TVector3 fBeamCenter; //Center of the beam with respect to the target center
  TVector3 fBeamPosition; //Position of the beam within the event (event-by-event)
  double fBeamAngularSpread; //FWHM of the angular spread of the beam
  double fBeamPositionSpread; //FWHM of the position spread of the beam
  double fBeamEnergySpread; //FWHM of the energy spread of the beam
  //
  
  //
  //Reaction Data
  std::string fPhysicsModelName;
  std::string fPhysicsConfigFileName;
  UNISIon fTheTarget;
  std::string fTargetMaterial;
  double fTargetThickness;
  double fTargetEffectiveThickness;
  double fTargetTilt;
  double fInteractionDistance; //distance from the entrance of the target at which the interaction occurs
  //
  
  //
  //Input Data
  TTree * fInputTree;
  TFile * fInputFile;
  TTreeReader * fTheReader;
  TTreeReaderValue<UNISRootEvent> * fTheInputEvent;
  //
  
  //
  //Output Data
  TTree * fTheTree;
  TFile * fTheFile;
  //
  
  //
  //Physics constants
  const double fUMAToMeV;
  const double fcmToum;
  //
  int fNumEvents; //Number of framework iterations
  //
  
  //
  //The detection setup
  UNISDetectionSetup * fExpSetup; 
  //
  
  //
  //The root event tools
  UNISRootEvent * fevt;
  //
  
  //
  //Modules
  UNISEventGenerator * fTheEventGenerator;
  //
  
  //
  //Methods used by the configuration process to read the input file
  int ProcessSetCommand(const char *);
  int ProcessDefineCommand(const char *);
  int ProcessAddCommand(const char *);
  //
  
  //
  //Calculation methods
  void GenerateBeam(); //Sets the fTheBeam UNISIon structure (called event-by-event)
  void InitializeTarget(); //Sets the fTheTarget UNISIon structure (called only once at the beginning of the simulation)
  void DetectEvent(); //Detect the event by using the detection setup from the UNISRootEvent data structure (after calling RegisterEvent())
  //
  
  //
  //Graphic tools
  TRint * fApp;
  void DisplayTracks(); //Display all the produced tracks as in fevt
  //
  
  //
  //Utility methods
  void PrintPercentage(int, int) const;
  //
};

#endif
