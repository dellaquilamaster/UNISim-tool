#ifndef UNISSEQUENTIALDECAYTWOBODY_H
#define UNISSEQUENTIALDECAYTWOBODY_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <UNISEventGenerator.h>
#include <UNISIon.h>
#include <UNISAngularDistribution.h>
#include <UNISShared.h>

class UNISSequentialDecayTwoBody : public UNISEventGenerator
{
public :
  UNISSequentialDecayTwoBody();
  ~UNISSequentialDecayTwoBody();
  
  int LoadConfiguration(const char *) override;
  
  std::vector<UNISIon> GetEvent() override;
  
private :
  TGenPhaseSpace * fRootGenerator;
  std::map <int, UNISIon *> fTheReactionProducts;
  std::map <int, UNISAngularDistribution *> fTheAngularDistributions;
  int ProcessSetCommand(const char *);
  int ProcessDefineCommand(const char *);
  int ReadAngularDistribution(const char * file_name, TH1D *);
  
  //
  // Distribution function used to produce excited states
  TF1 * fBreitWignerDistribution;
  //
  
  //
  //Method to simulate the decay of a secondary particle (recursive)
  void SecondaryDecay(UNISIon *, std::vector<UNISIon> &, UNISAngularDistribution * ang_distr=0);
  //
  
};

#endif
