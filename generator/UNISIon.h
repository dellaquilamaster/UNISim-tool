#ifndef UNISION_H
#define UNISION_H

#include <map>
#include <TLorentzVector.h>

typedef struct ExcitedState_
{
  double Ex;
  double Gamma;
} ExcitedState;

class UNISIon
{
public :
  TLorentzVector fMomentum;
  int fZ; //Charge
  int fA; //Mass Number
  double fMass; //Mass MeV/c2 (does not account for the excitation energy)
  std::vector<ExcitedState> fSpectroscopy; //Excited states
  std::map<int, UNISIon *> fSecondaryParticles; //Secondary particles generated by its decay
};


#endif
