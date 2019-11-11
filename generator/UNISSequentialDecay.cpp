#include <UNISSequentialDecay.h>

//____________________________________________________
UNISSequentialDecay::UNISSequentialDecay() :
fRootGenerator(new TGenPhaseSpace()),
fBreitWignerDistribution(new TF1("BreitWignerDistribution", "1/((x-[0])^2+([1]^2)/4)", -100, 100))
{
  fBreitWignerDistribution->SetNpx(500);
}

//____________________________________________________
UNISSequentialDecay::~UNISSequentialDecay()
{
  delete fRootGenerator;
  delete fBreitWignerDistribution;
}


//____________________________________________________
int UNISSequentialDecay::LoadConfiguration(const char * file_name)
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
    
  return NRead;
}

//____________________________________________________
std::vector<UNISIon> UNISSequentialDecay::GetEvent()
{
  std::vector<UNISIon> AnEvent;
  
  //
  //Generating the main particles emitted in the first step of the reaction
  //Then, the decay of each of those 2 particles will be simulated individually, by following eventually the whole decay chain
  //of each of the emitted fragments
  //
  //Calculating the total momentum available in the reaction
  TLorentzVector TotalMomentum = fTheBeam.fMomentum+fTheTarget.fMomentum;
  //
  const int NumParticles = fTheReactionProducts.size(); //number of particles in the primary reaction
  double Masses[NumParticles]; //Masses of each fragment in the exit channel
  //
  //Initializing masses
  for(int i=0; i<NumParticles; i++) {
    Masses[i]=fTheReactionProducts[i]->fMass;
    //
    //Choose the excited state to populate in the i-th fragment (if any)
    if(fTheReactionProducts[i]->fSpectroscopy.size()>0) {
      const int WhichState = int(fRandomGen->Uniform(0,fTheReactionProducts[i]->fSpectroscopy.size()));
      double ExcitationEnergy;
      if(fTheReactionProducts[i]->fSpectroscopy[WhichState].Gamma>0) {
        fBreitWignerDistribution->SetParameters(fTheReactionProducts[i]->fSpectroscopy[WhichState].Ex, fTheReactionProducts[i]->fSpectroscopy[WhichState].Gamma);
        ExcitationEnergy = fBreitWignerDistribution->GetRandom();
      } else {
        ExcitationEnergy = fTheReactionProducts[i]->fSpectroscopy[WhichState].Ex;
      }
      Masses[i]+=ExcitationEnergy;      
    }
    //
    
    //NOTE
//     printf("Particella primaria %d -> m=%f\n", i, Masses[i]);
    //NOTE
    
  }
  //
  //Simulation of the primary reaction
  fRootGenerator->SetDecay(TotalMomentum,NumParticles,Masses);
  fRootGenerator->Generate();
  //
  //Setting momenta of each particle
  for(int i=0; i<NumParticles; i++) {
    TLorentzVector * TheDecayFragment = fRootGenerator->GetDecay(i); //Getting the i-th particle in the decay
    fTheReactionProducts[i]->fMomentum=*TheDecayFragment;
  }
  //
  
  //
  //Calculate (if needed) the secondary decay of each particle in the final state
  //And write the result on the global event
  //This method is recursive as it will calculate the secondary decay of 
  //each of the secondary particles and so on
  for(int i=0; i<NumParticles; i++)
  {
    if(fTheReactionProducts[i]->fSecondaryParticles.size()==0) {
      AnEvent.push_back(*fTheReactionProducts[i]); //The particle does not decay, writing particle in the event
    } else SecondaryDecay(fTheReactionProducts[i],AnEvent); //If the particle has a decay pattern, process the decay pattern recursively till the end of the chain
  }
  //
  
  return AnEvent;
}

//The method SecondaryDecay takes care of calculating the secondary decay of 
//all particles originating by an initial particle in the exit channel and
//write all of them in the event used as the output of the GetEvent method.
//____________________________________________________
void UNISSequentialDecay::SecondaryDecay(UNISIon * TheParticle, std::vector<UNISIon> & TheEvent)
{
  //Setting the decay of the current step
  const int NumSecondaryParticles = TheParticle->fSecondaryParticles.size();
  double Masses[NumSecondaryParticles];
  for(int i=0; i<NumSecondaryParticles; i++) {
    Masses[i]=(TheParticle->fSecondaryParticles)[i]->fMass;
    //
    //Choose the excited state to populate in the i-th fragment (if any)
    if((TheParticle->fSecondaryParticles)[i]->fSpectroscopy.size()>0) {
      const int WhichState = int(fRandomGen->Uniform(0,(TheParticle->fSecondaryParticles)[i]->fSpectroscopy.size()));
      double ExcitationEnergy;
      if((TheParticle->fSecondaryParticles)[i]->fSpectroscopy[WhichState].Gamma>0) {
        fBreitWignerDistribution->SetParameters((TheParticle->fSecondaryParticles)[i]->fSpectroscopy[WhichState].Ex, (TheParticle->fSecondaryParticles)[i]->fSpectroscopy[WhichState].Gamma);
        ExcitationEnergy = fBreitWignerDistribution->GetRandom();
      } else {
        ExcitationEnergy = (TheParticle->fSecondaryParticles)[i]->fSpectroscopy[WhichState].Ex; 
      }
      Masses[i]+=ExcitationEnergy;
    }
    //
    
    //NOTE
//     printf("Particella secondaria %d derivante dalla particella di massa %f -> m=%f\n", i, TheParticle->fMass, Masses[i]);
    //NOTE

  }
  //
  //Setting the generator of the secondary decay
  fRootGenerator->SetDecay(TheParticle->fMomentum,NumSecondaryParticles,Masses);
  //
  //Processing the secondary decay
  fRootGenerator->Generate();
  //
  //Setting momenta of each particle
  for(int i=0; i<NumSecondaryParticles; i++) {
    TLorentzVector * TheDecayFragment = fRootGenerator->GetDecay(i); //Getting the i-th particle in the decay
    (TheParticle->fSecondaryParticles)[i]->fMomentum=*TheDecayFragment; 
  }
  //
  
  //
  //Writing particles to the event or process additional secondary decays (recursively)
  for(int i=0; i<NumSecondaryParticles; i++)
  {
    if((TheParticle->fSecondaryParticles)[i]->fSecondaryParticles.size()==0) {
      TheEvent.push_back(*(TheParticle->fSecondaryParticles)[i]); //The particle does not decay, writing particle in the event
    } else SecondaryDecay((TheParticle->fSecondaryParticles)[i],TheEvent); //If the particle has a decay pattern, process the decay pattern recursively till the end of the chain
  }
  //
  
  return;
}

//____________________________________________________
int UNISSequentialDecay::ProcessSetCommand(const char * line)
{
  std::string InputLine(line);
  std::istringstream LineStream(InputLine);
  
  std::string Command;
  std::string WhatToSet;
  
  LineStream >> Command >> WhatToSet;
  
  if(WhatToSet.compare("spectroscopy")==0) {
    //Setting a new excited state to a certain nucleus
    std::string NucleusToSet;
    LineStream>>NucleusToSet;
    std::istringstream TheStream(NucleusToSet);
    std::string ANucleus;
    std::vector<int> PathToTheNucleus;
    //NucleusToSet contains the number of particle and sub-particle to which to set the excited states
    //We need to separate it by using "_" as the delimiter character
    while(std::getline(TheStream, ANucleus, '_')) {
      PathToTheNucleus.push_back(std::stoi(ANucleus)); 
    }
    //
    //Retrieving the quantities to set
    std::string TheSpectroscopyQuantity;
    ExcitedState NewState;
    while(LineStream>>TheSpectroscopyQuantity) {
      if(TheSpectroscopyQuantity.find("-Ex=")!=std::string::npos) {
        NewState.Ex=std::stof(TheSpectroscopyQuantity.substr(TheSpectroscopyQuantity.find("-Ex=")+4));
      } else if(TheSpectroscopyQuantity.find("-Gamma=")!=std::string::npos) {
        NewState.Gamma=std::stof(TheSpectroscopyQuantity.substr(TheSpectroscopyQuantity.find("-Gamma=")+7));
      }
    }
    //
    //Setting quantities to the correct nucleus
    UNISIon * TheNucleusToSet=fTheReactionProducts[PathToTheNucleus[0]];
    for(int i=1; i<PathToTheNucleus.size(); i++) {
      TheNucleusToSet=TheNucleusToSet->fSecondaryParticles[PathToTheNucleus[i]];
    }
    
    TheNucleusToSet->fSpectroscopy.push_back(NewState);
    //    
  } else if(WhatToSet.compare("decay")==0) {
    //Setting a new decay fragment of a certain nucleus
    std::string NucleusToSet;
    LineStream>>NucleusToSet;
    std::istringstream TheStream(NucleusToSet);
    std::string ANucleus;
    std::vector<int> PathToTheNucleus;
    //NucleusToSet contains the number of particle and sub-particle to which to set decay
    //We need to separate it by using "_" as the delimiter character
    while(std::getline(TheStream, ANucleus, '_')) {
      PathToTheNucleus.push_back(std::stoi(ANucleus)); 
    }
    //
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
    //
    //Setting quantities to the correct nucleus
    UNISIon * TheNucleusToSet=fTheReactionProducts[PathToTheNucleus[0]];
    for(int i=1; i<PathToTheNucleus.size(); i++) {
      TheNucleusToSet=TheNucleusToSet->fSecondaryParticles[PathToTheNucleus[i]];
    }
    //        
    TheNucleusToSet->fSecondaryParticles[ParticleToSet]=NewParticle;
  } else {
    return 0; 
  }
  
  return 1;
}

//____________________________________________________
int UNISSequentialDecay::ProcessDefineCommand(const char * line)
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
