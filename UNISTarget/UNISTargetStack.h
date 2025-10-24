#ifndef UNISTARGETSTACK_H
#define UNISTARGETSTACK_H

//Not having active target at all should be handled by framework

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <UNISShared.h>
#include <UNISIon.h>
#include <UNISTarget.h>

#include <EnergyLossModule.h>
#include <nuclear_masses.h>

typedef struct Interaction_Point_
{
  unsigned int target_number;
  double depth;
} Interaction_Point;

class UNISTargetStack{
public :
  UNISTargetStack();
  ~UNISTargetStack();
  
  void SetName(const char *);

  void ListStack();
  void AddTarget(UNISTarget* target);
  void RemoveTarget(UNISTarget* target);
  UNISTarget GetTarget_fromStack(unsigned int target_number=0);
  void SetStackTiltX(double rotationAround_xAxis); //in degrees
  void SetStackTiltY(double rotationAround_yAxis); //in degrees

  int GetNum_ofTargets_inStack();
  const TVector3 GetNormalVector_ofStack() const;
  Interaction_Point GetInteractionPoint();
  //If the beam stopps in material, returned UNISIon with fMomentum (0, 0, 0, mass_gs)
  UNISIon GetBeam_preInteraction(const UNISIon& preTargetBeam, Interaction_Point point_of_interaction); 
  //If the particle stopps in material, returned UNISIon with fMomentum (0, 0, 0, mass_gs)
  std::vector<UNISIon> PropagateParticles(const std::vector<UNISIon>& Particles, Interaction_Point point_of_interaction);

protected :
  TRandom3 * fRandomGen;
  std::string fName;
  double total_thickness;
  std::vector<unsigned int> active_targets_indices;
  std::vector<double> active_targets_limits;
  //std::map<std::pair<unsigned int, unsigned int>, double> active_targets; //index (regarding all active targets), num of target (regarding all targets), boundary 
  std::vector<UNISTarget> target_stack;
  TVector3 norm_vector; //vector perpendicular to the target surface
  TVector3 tan_vector; //vector parallel to the target surface
};

#endif
