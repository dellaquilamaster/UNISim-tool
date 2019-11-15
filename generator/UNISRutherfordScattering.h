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
#include <UNISShared.h>

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
  double fMinAngle; //minimum scattering angle used to generate elastic scattering data
  double fMaxAngle; //maximum scattering angle used to generate elastic scattering data
  TF1 * fAngularDistribution;
  bool fIsDistributionUniform; //true if the user chooses an uniform distribution instead of a rutherford-like
  //
};

#endif
