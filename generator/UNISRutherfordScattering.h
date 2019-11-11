#ifndef UNISRUTHERFORDSCATTERING_H
#define UNISRUTHERFORDSCATTERING_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>
#include <TMath.h>
#include <TF1.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <UNISEventGenerator.h>
#include <UNISIon.h>
#include <RelativisticKinematics.h>
#include <shared.h>

class UNISRutherfordScattering : public UNISEventGenerator
{
public :
  UNISRutherfordScattering();
  ~UNISRutherfordScattering();
  
  int LoadConfiguration(const char *) override;
  
  std::vector<UNISIon> GetEvent() override;
  
private :
  RelativisticKinematics * fTheTwoBodyKinematicsModule;
  std::map <int, UNISIon *> fTheReactionProducts;
  int ProcessDefineCommand(const char *);
  int ProcessSetCommand(const char *);
  
  //
  // Distribution function used to generate the scattering angle
  double fMinAngle; //minimum scattering angle used to generate Rutherford scattering data
  TF1 * fRutherfordDistribution;
  //
};

#endif
