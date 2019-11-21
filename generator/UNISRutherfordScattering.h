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
  bool fIsUserDefinedDistribution; //true if the user chooses a customized angular distribution instead of a rutherford-like
  std::string fUserDefinedDistrubutionFormula; //the formula for the user-defined angular distribution
  int fUserDefinedDistributionNumParameters; //number of parameters of the user-defined angular distribution
  double * fUserDefinedDistributionParameters; //the parameters of the user-defined angular distribution
  //
};

#endif
