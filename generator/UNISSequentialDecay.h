#ifndef UNISSEQUENTIALDECAY_H
#define UNISSEQUENTIALDECAY_H

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
#include <UNISShared.h>

class UNISSequentialDecay : public UNISEventGenerator
{
public :
  UNISSequentialDecay();
  ~UNISSequentialDecay();
  
  int LoadConfiguration(const char *) override;
  
  std::vector<UNISIon> GetEvent() override;
  
private :
  TGenPhaseSpace * fRootGenerator;
  std::map <int, UNISIon *> fTheReactionProducts;
  int ProcessSetCommand(const char *);
  int ProcessDefineCommand(const char *);
  
  //
  // Distribution function used to produce excited states
  TF1 * fBreitWignerDistribution;
  //
  
  //
  //Method to simulate the decay of a secondary particle (recursive)
  void SecondaryDecay(UNISIon *, std::vector<UNISIon> &);
  //
};

#endif
