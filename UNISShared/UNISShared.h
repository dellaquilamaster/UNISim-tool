#ifndef SHARED_H
#define SHARED_H

class EnergyLossModule;
class nuclear_masses;

//Maximum multiplicity of an event
extern const int gMaxEvtMulti;
//Nuclear masses database
extern nuclear_masses * gNucData;
//Module for Energy Loss
extern EnergyLossModule * gLISEELossModule;
//Random number seed
extern double gRandomSeed;

#endif
