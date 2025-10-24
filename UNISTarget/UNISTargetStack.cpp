#include "UNISTargetStack.h"

//________________________________________________
UNISTargetStack::UNISTargetStack() :
fRandomGen(new TRandom3 (gRandomSeed)),
fName(""),
total_thickness(0),
norm_vector(0, 0, 1),
tan_vector(1, 0, 0)
{}

//________________________________________________
UNISTargetStack::~UNISTargetStack()
{
  delete fRandomGen;
}

//________________________________________________
void UNISTargetStack::SetName(const char * target_stack_name)
{
  fName.assign(target_stack_name);
}

//_______________________________________
void UNISTargetStack::AddTarget(UNISTarget* target) 
{  
  if (!target) return;

  target_stack.push_back(*target);
  
  if(target->TargetIsActive())
  {
    total_thickness += target->GetTargetThickness_um();
    active_targets_indices.push_back(target_stack.size() - 1);
    active_targets_limits.push_back(total_thickness);
    //active_targets[{active_targets.size(), target_stack.size() - 1}] = total_thickness;
  }
}

//_______________________________________
void UNISTargetStack::RemoveTarget(UNISTarget* target) 
{
  if (!target) return;

  std::vector<unsigned int> active_targets_indices_new;
  unsigned int total_index = 0;  // counts all targets

  // If target inactive
  if (!target->TargetIsActive()) {
    for (auto it = target_stack.begin(); it != target_stack.end(); ) 
    {
      if (*it == *target) 
      {
        target_stack.erase(it);
        continue; // only one target to remove
      }
      if(it->TargetIsActive())
      {
        active_targets_indices_new.push_back(total_index);
      }
      ++it;
      ++total_index;
    }
    active_targets_indices = std::move(active_targets_indices_new);
    return;
  }

  // If target is active
  double total_thickness_new = 0.0;
  std::vector<double> active_targets_limits_new;

  //std::map<std::pair<unsigned int, unsigned int>, double> active_targets_new;

  //unsigned int active_index = 0; // counts only active targets

  for (auto it = target_stack.begin(); it != target_stack.end(); ) {
    if (*it == *target)
    {
      it = target_stack.erase(it); // safe erase
      continue;
    } 

    if (it->TargetIsActive()) 
    {
      total_thickness_new += it->GetTargetThickness_um();
      active_targets_indices_new.push_back(total_index);
      active_targets_limits_new.push_back(total_thickness_new);
      //active_targets_new[{active_index, total_index}] = total_thickness_new;
      //++active_index;
    }

    ++total_index;
    ++it;
  }

  total_thickness = total_thickness_new;
  active_targets_limits = std::move(active_targets_limits_new);
  active_targets_indices = std::move(active_targets_indices_new);
  //active_targets = std::move(active_targets_new); //faster than "="; doesn't copy each element
}

//_______________________________________
UNISTarget UNISTargetStack::GetTarget_fromStack(unsigned int target_number) {
  if (target_number >= target_stack.size()) return UNISTarget("", "", false, 0, 0);

  return target_stack[target_number];
}

//_______________________________________
void UNISTargetStack::ListStack() 
{
  std::cout << "===== Target Stack: " << fName << " =====" << std::endl;
  for (size_t i = 0; i < target_stack.size(); ++i) {
    std::cout << "Layer " << i << ":" << std::endl;
    target_stack[i].Print();
    std::cout << std::endl;
  }
  std::cout << "=========================================" << std::endl;
}

//_______________________________________
void UNISTargetStack::SetStackTiltX(double rotationAround_xAxis) 
{
  norm_vector.RotateX(TMath::DegToRad()*rotationAround_xAxis);
  tan_vector.RotateX(TMath::DegToRad()*rotationAround_xAxis);
}

//_______________________________________
void UNISTargetStack::SetStackTiltY(double rotationAround_yAxis) 
{
  norm_vector.RotateY(TMath::DegToRad()*rotationAround_yAxis);
  tan_vector.RotateY(TMath::DegToRad()*rotationAround_yAxis);
}

//_______________________________________
const TVector3 UNISTargetStack::GetNormalVector_ofStack() const 
{
  return norm_vector;
}

//_______________________________________
int UNISTargetStack::GetNum_ofTargets_inStack()
{
  return target_stack.size();
}


//_______________________________________
Interaction_Point UNISTargetStack::GetInteractionPoint()
{
  ////choose randomly from active targets, but weighted according to their thicknesses
  //std::vector<double> target_bounds;
  //std::vector<int> target_indices;
  //double total = 0.0;

  //for (size_t i = 0; i < target_stack.size(); ++i)
  //{
  //  if(!target_stack[i].TargetIsActive()) continue; //if target is not active, continue

  //  total += target_stack[i].GetTargetThickness_um();
  //  target_bounds.push_back(total);
  //  target_indices.push_back(i);
  //}

  if(target_stack.size() == 0)
    return Interaction_Point{0, 0.0};

  double rand_point = fRandomGen->Uniform(0.0, total_thickness); //returns 0 if thickess is 0
  unsigned int target_num = 0;
  double depth_inTarget = 0.0;
  double prev_boundary;

  if(total_thickness != 0)
  {
    //for(const auto& one_active_target : active_targets) //active_targets[{activ_index, tot_index}] = boundary
    for(size_t i = 0; i<active_targets_limits.size(); ++i)
    {
      if(rand_point < active_targets_limits[i])
      {
        target_num = active_targets_indices[i];
        prev_boundary = (i == 0) ? 0.0 : active_targets_limits[i-1];
        depth_inTarget = rand_point - prev_boundary;
        break;
      }
    }
  }
  else
  {
    unsigned int rand_target_num = fRandomGen -> Integer(active_targets_limits.size());
    target_num = active_targets_indices[rand_target_num];
  }
  
  Interaction_Point point;
  point.target_number = target_num;
  point.depth = depth_inTarget;

  return point;
}

//_______________________________________
UNISIon UNISTargetStack::GetBeam_preInteraction(const UNISIon& preTargetBeam, Interaction_Point point_of_interaction)
{
  UNISIon propagatedBeam = preTargetBeam;
  unsigned int reactionTarget = point_of_interaction.target_number;
  double reactionDepth = point_of_interaction.depth;

  int beamZ = propagatedBeam.fZ;
  int beamA = propagatedBeam.fA;
  double beam_mass = propagatedBeam.fMass;
  TLorentzVector beamV = propagatedBeam.fMomentum;
  
  TVector3 beam_momentum = propagatedBeam.fMomentum.Vect(); // gets the 3-momentum
  TVector3 beam_direction = beam_momentum.Mag() > 0 ? beam_momentum.Unit() : TVector3(0, 0, 0);
  TVector3 stack_direction = norm_vector.Mag() > 0 ? norm_vector.Unit() : TVector3(0, 0, 0);
  
  double kinE_preReaction = beamV.E() - beam_mass;
  double eLoss;

  std::string material;
  double thickness_um, eff_thickness, dotProduct;

  if( reactionTarget == 0 && reactionDepth == 0)
    return propagatedBeam;

  //Energy loss in targets up to the target where reaction happend
  for (size_t i = 0; i < reactionTarget; ++i)
  {
    material = target_stack[i].GetTargetMaterial();
    thickness_um = target_stack[i].GetTargetThickness_um();
    if(thickness_um == 0) continue;
    dotProduct = beam_direction.Dot(stack_direction);
    eff_thickness = dotProduct > 1e-12 ? thickness_um/(dotProduct) : 9999.9;
    eLoss = gLISEELossModule->GetEnergyLoss(beamZ, beamA, kinE_preReaction, material.c_str(), eff_thickness);
    kinE_preReaction = kinE_preReaction - eLoss;
  }

  //Energy loss in target, up to the point of reaction
  material = target_stack[reactionTarget].GetTargetMaterial();
  thickness_um = target_stack[reactionTarget].GetTargetThickness_um();
  //debug
  //if(reactionDepth > thickness_um)
  //{
  //  std::cerr << "[GetBeamEnergy_preInteraction] Reaction Depth is bigger than Thickness of Target in um!" << std::endl;
  //  propagatedBeam.fMomentum = TLorentzVector(0, 0, 0, -1);
  //  return propagatedBeam;    
  //}

  if(reactionDepth == 0)
    kinE_preReaction = kinE_preReaction;
  else
  {
    dotProduct = beam_direction.Dot(stack_direction);
    eff_thickness = dotProduct > 1e-12 ? (reactionDepth)/(dotProduct) : 9999.9;
    eLoss = gLISEELossModule->GetEnergyLoss(beamZ, beamA, kinE_preReaction, material.c_str(), eff_thickness);
    kinE_preReaction = kinE_preReaction - eLoss;    
  }

  double total_energy = kinE_preReaction + beam_mass;
  double momentum = sqrt(pow(total_energy, 2) - pow(beam_mass, 2));
  double theta = beam_direction.Theta();
  double phi = beam_direction.Phi();
  TLorentzVector q(momentum*sin(theta)*cos(phi),momentum*sin(theta)*sin(phi),momentum*cos(theta),total_energy);

  propagatedBeam.fMomentum = q;

  return propagatedBeam;
}

//_______________________________________
std::vector<UNISIon> UNISTargetStack::PropagateParticles(const std::vector<UNISIon>& Particles, Interaction_Point point_of_interaction)
{
  std::vector<UNISIon> propagatedParticles = Particles;
  unsigned int reactionTarget = point_of_interaction.target_number;
  double reactionDepth = point_of_interaction.depth;

  std::string material;
  int Z, A;
  double thickness_um, eff_thickness, eLoss, kinE, total_energy, momentum, theta, phi, dotProduct;
  TLorentzVector q_LV;
  TVector3 p_mom, p_direction;
  TLorentzVector q_final;

  TVector3 stack_direction = norm_vector.Mag() > 0 ? norm_vector.Unit() : TVector3(0, 0, 0);

  for(auto& particle : propagatedParticles)
  {
    //Extracting info from particles
    Z = particle.fZ;
    A = particle.fA;
    //mass = particle.fMass;
    q_LV = particle.fMomentum;
    p_mom = particle.fMomentum.Vect(); // gets the 3-momentum
    p_direction = p_mom.Mag() > 0 ? p_mom.Unit() : TVector3(0, 0, 0);
    kinE = q_LV.E() - q_LV.M();

    //Differentiate if particle goes to theta < 90 or theta > 90
    if(p_direction.Theta() < (TMath::Pi()/2.0)) // goes "foreward"
    {
      material = target_stack[reactionTarget].GetTargetMaterial();
      thickness_um = target_stack[reactionTarget].GetTargetThickness_um();

      if(thickness_um == 0) kinE = kinE;
      else
      {
        dotProduct = p_direction.Dot(stack_direction);
        eff_thickness = fabs(dotProduct) > 1e-12 ? (thickness_um - reactionDepth)/fabs(dotProduct) : 9999.9;
        eLoss = gLISEELossModule->GetEnergyLoss(Z, A, kinE, material.c_str(), eff_thickness, 4);
        kinE = kinE - eLoss;        
      }

      //Energy loss in target, from the point of reaction
      //#debug
      //if(reactionDepth > thickness_um)
      //{
      //  std::cerr << "[PropagateParticles] Reaction Depth is bigger than Thickness of Target in um!" << std::endl;
      //}

      if(reactionTarget < (target_stack.size() - 1))
      {
        for (size_t i = reactionTarget + 1 ; i < target_stack.size(); ++i)
        {
          material = target_stack[i].GetTargetMaterial();
          thickness_um = target_stack[i].GetTargetThickness_um();
          if(thickness_um == 0) continue;
          dotProduct = p_direction.Dot(stack_direction);
          eff_thickness = fabs(dotProduct) > 1e-12 ? thickness_um/fabs(dotProduct) : 9999.9;
          eLoss = gLISEELossModule->GetEnergyLoss(Z, A, kinE, material.c_str(), eff_thickness, 4);
          kinE = kinE - eLoss;
        }
      }
    }
    else
    {
      material = target_stack[reactionTarget].GetTargetMaterial();
      thickness_um = target_stack[reactionTarget].GetTargetThickness_um();

      //Energy loss in target, from the point of reaction
      //if(reactionDepth > thickness_um)
      //{
      //  std::cerr << "[PropagateParticles] Reaction Depth is bigger than Thickness of Target in um!" << std::endl;
      //}

      if(reactionDepth == 0) kinE = kinE;
      else
      {
        dotProduct = p_direction.Dot(stack_direction);
        eff_thickness = fabs(dotProduct) > 1e-12 ? (reactionDepth)/fabs(dotProduct) : 9999.9;
        eLoss = gLISEELossModule->GetEnergyLoss(Z, A, kinE, material.c_str(), eff_thickness, 4);
        kinE = kinE - eLoss;        
      }

      if(reactionTarget > 0)
      {
        for (int i = static_cast<int>(reactionTarget) - 1; i >= 0; --i)
        {
          material = target_stack[i].GetTargetMaterial();
          thickness_um = target_stack[i].GetTargetThickness_um();
          if(thickness_um == 0) continue;
          dotProduct = p_direction.Dot(stack_direction);
          eff_thickness = fabs(dotProduct) > 1e-12 ? thickness_um/fabs(dotProduct) : 9999.9;
          eLoss = gLISEELossModule->GetEnergyLoss(Z, A, kinE, material.c_str(), eff_thickness, 4);
          kinE = kinE - eLoss;
        }
      }
    }

    total_energy = kinE + q_LV.M();
    momentum = sqrt(pow(total_energy, 2) - pow(q_LV.M(), 2));
    theta = p_direction.Theta();
    phi = p_direction.Phi();
    TLorentzVector q_final(momentum*sin(theta)*cos(phi),momentum*sin(theta)*sin(phi),momentum*cos(theta),total_energy);

    particle.fMomentum = q_final;
  }

  return propagatedParticles;
}