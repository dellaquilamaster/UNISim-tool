#ifndef SHARED_H
#define SHARED_H

#include <nuclear_masses.h>
#include <EnergyLossModule.h>

//Maximum multiplicity of an event
extern const int gMaxEvtMulti;
//Nuclear masses database
extern nuclear_masses * gNucData;
//Module for Energy Loss
extern EnergyLossModule * gLISEELossModule;
//Random number seed
extern double gRandomSeed;

#endif
