#ifndef UNISEVENTGENERATOR_H
#define UNISEVENTGENERATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <nuclear_masses.h>
#include <UNISShared.h>
#include <UNISIon.h>

class UNISEventGenerator
{
public :
  UNISEventGenerator();
  ~UNISEventGenerator();
  
  virtual int LoadConfiguration(const char *);
  
  void SetBeam(UNISIon &);
  void SetTarget(UNISIon &);
  
  virtual std::vector<UNISIon> GetEvent();
  
protected :
  TRandom3 * fRandomGen;
  UNISIon fTheBeam;
  UNISIon fTheTarget;
  
};

#endif
