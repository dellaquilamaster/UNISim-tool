#include <UNISRutherfordScattering.h>

//____________________________________________________
UNISRutherfordScattering::UNISRutherfordScattering() :
fTheTwoBodyKinematicsModule(new RelativisticKinematics()),
fMinAngle(0)
{}

//____________________________________________________
UNISRutherfordScattering::~UNISRutherfordScattering()
{
  delete fTheTwoBodyKinematicsModule;
  delete fRutherfordDistribution;
}


//____________________________________________________
int UNISRutherfordScattering::LoadConfiguration(const char * file_name)
{
  std::ifstream FileIn(file_name);
  if(!FileIn.is_open()) {
    return -1;
  }
  
  int NRead=0;

  while (!FileIn.eof())
  {
    std::string LineRead;
    std::getline(FileIn, LineRead);

    std::string LineReadCommentLess (LineRead.substr(0, LineRead.find("*")));

    if(LineReadCommentLess.empty()) continue;
    if(LineReadCommentLess.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineReadCommentLess);
    
    std::string Command;
    
    LineStream >> Command;
    
    if(Command.compare("set")==0) {
      NRead+=ProcessSetCommand(LineReadCommentLess.c_str());
    } else if (Command.compare("define")==0) {
      NRead+=ProcessDefineCommand(LineReadCommentLess.c_str());
    }
  }
  
  //
  //Configuring the angular distribution
  fRutherfordDistribution = new TF1("RutherfordDistribution", "(([0]*[1]*(1./137)*197)/(4*[2]*(sin(x/2.))^2))^2", fMinAngle, TMath::Pi());
  fRutherfordDistribution->SetParameter(0, fTheReactionProducts[0]->fZ); //Charge of the projectile (ejectile)
  fRutherfordDistribution->SetParameter(1, fTheReactionProducts[1]->fZ); //Charge of the target (recoil)
  fRutherfordDistribution->SetNpx(500);
  //
  
  return NRead;
}

//____________________________________________________
std::vector<UNISIon> UNISRutherfordScattering::GetEvent()
{
  std::vector<UNISIon> AnEvent;
  
  //
  //Generating the main particles emitted in the first step of the reaction
  //Then, the decay of each of those 2 particles will be simulated individually, by following eventually the whole decay chain
  //of each of the emitted fragments
  //
  //Kinetic energy of the beam
  double BeamEnergy = fTheBeam.fMomentum.E()-fTheBeam.fMomentum.M();
  //Calculating the total momentum available in the reaction
  TLorentzVector TotalMomentum = fTheBeam.fMomentum+fTheTarget.fMomentum;
  //Setting the energy to the Rutherford Distribution
  fRutherfordDistribution->SetParameter(2, BeamEnergy); //kinetic energy of the beam
  //
  const int NumParticles = fTheReactionProducts.size(); //number of particles in the primary reaction
  double Masses[NumParticles]; //Masses of each fragment in the exit channel
  //
  //Initializing masses
  for(int i=0; i<NumParticles; i++) {
    Masses[i]=fTheReactionProducts[i]->fMass;
  }
  //    
  //
  //Simulation of the primary reaction (elastic scattering according to the Rutherford Distribution)
  double ThetaRutherford = fRutherfordDistribution->GetRandom();
  double PhiEjectile = fRandomGen->Uniform(0,2*TMath::Pi());
  double EjectileEnergy = fTheTwoBodyKinematicsModule->GetEnergyTwoBodyEjectile(ThetaRutherford, BeamEnergy/gNucData->get_mass_Z_A_uma(fTheBeam.fZ,fTheBeam.fA), gNucData->get_mass_Z_A_uma(fTheBeam.fZ,fTheBeam.fA), gNucData->get_mass_Z_A_uma(fTheTarget.fZ,fTheTarget.fA));
  double EjectileTotalEnergy = EjectileEnergy+fTheBeam.fMass;
  double EjectileMomentumModule = sqrt(pow(EjectileTotalEnergy,2)-pow(fTheBeam.fMass,2));
  TLorentzVector TheEjectile (EjectileMomentumModule*sin(ThetaRutherford)*cos(PhiEjectile),EjectileMomentumModule*sin(ThetaRutherford)*sin(PhiEjectile),EjectileMomentumModule*cos(ThetaRutherford),EjectileTotalEnergy);
  TLorentzVector TheRecoil = TotalMomentum - TheEjectile;
  
  //
  //Setting momenta of each particle
  for(int i=0; i<NumParticles; i++) {
    if(fTheReactionProducts[i]->fZ==fTheBeam.fZ && fTheReactionProducts[i]->fA==fTheBeam.fA) {
      fTheReactionProducts[i]->fMomentum=TheEjectile;
      AnEvent.push_back(*fTheReactionProducts[i]);
    } else if (fTheReactionProducts[i]->fZ==fTheTarget.fZ && fTheReactionProducts[i]->fA==fTheTarget.fA) {
      fTheReactionProducts[i]->fMomentum=TheRecoil;
      AnEvent.push_back(*fTheReactionProducts[i]);
    }
  }
  //
  
  return AnEvent;
}

//____________________________________________________
int UNISRutherfordScattering::ProcessDefineCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);  
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream>>Command>>WhatToSet;
  
  if(WhatToSet.compare("particle")==0) {
    int ParticleToSet;
    int Charge;
    int Mass;
    LineStream>>ParticleToSet;
    while(LineStream>>WhatToSet) {
      if(WhatToSet.find("-Z=")!=std::string::npos) {
        Charge=std::stoi(WhatToSet.substr(WhatToSet.find("-Z=")+3));
      } else if(WhatToSet.find("-A=")!=std::string::npos) {
        Mass=std::stoi(WhatToSet.substr(WhatToSet.find("-A=")+3));
      }
    }
    //
    //Adding a new particle
    UNISIon * NewParticle = new UNISIon();
    NewParticle->fZ=Charge;
    NewParticle->fA=Mass;
    NewParticle->fMass=gNucData->get_mass_Z_A(Charge,Mass);
    
    fTheReactionProducts[ParticleToSet]=NewParticle;
    //
  } else {
    return 0; 
  }

  return 1;
}

//____________________________________________________
int UNISRutherfordScattering::ProcessSetCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);  
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream>>Command>>WhatToSet;
  
  if(WhatToSet.compare("min_angle")==0) {
    double TheAngle;
    LineStream>>TheAngle;
    fMinAngle=TheAngle*TMath::DegToRad();
  } else {
    return 0; 
  }

  return 1;
}

