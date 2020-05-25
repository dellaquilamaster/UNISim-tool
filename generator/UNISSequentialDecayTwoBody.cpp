#include <UNISSequentialDecayTwoBody.h>

//____________________________________________________
UNISSequentialDecayTwoBody::UNISSequentialDecayTwoBody() :
fRootGenerator(new TGenPhaseSpace()),
fBreitWignerDistribution(new TF1("BreitWignerDistribution", "1/((x-[0])^2+([1]^2)/4)", -100, 100))
{
  fBreitWignerDistribution->SetNpx(500);
}

//____________________________________________________
UNISSequentialDecayTwoBody::~UNISSequentialDecayTwoBody()
{
  delete fRootGenerator;
  delete fBreitWignerDistribution;
}


//____________________________________________________
int UNISSequentialDecayTwoBody::LoadConfiguration(const char * file_name)
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
std::vector<UNISIon> UNISSequentialDecayTwoBody::GetEvent()
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
  //Simulation of the primary reaction (2-body)
  TH1D * PrimaryDecayAngDistr = fTheAngularDistributions[0]->fAngularDistribution;
  double theta_cms = PrimaryDecayAngDistr ? PrimaryDecayAngDistr->GetRandom()*TMath::DegToRad() : acos(gRandom->Uniform(0,1));
  double phi = gRandom->Uniform(0,2*TMath::Pi());
  const double InvariantMass=TotalMomentum.M();
  const double MomentumCmModule=sqrt((pow(InvariantMass,4)+pow((pow(Masses[0],2)-pow(Masses[1],2)),2)-2*pow(InvariantMass,2)*(pow(Masses[0],2)+pow(Masses[1],2)))/(4*pow(InvariantMass,2))); //Module of the momentum in the cm
  TVector3 EjectileMomentumCm(MomentumCmModule*sin(theta_cms)*cos(phi),MomentumCmModule*sin(theta_cms)*sin(phi),MomentumCmModule*cos(theta_cms)); //Momentum of the ejectile in the cm
  TLorentzVector EjectileMomentumEnergy(EjectileMomentumCm,sqrt(pow(MomentumCmModule,2)+pow(Masses[0],2))); //Quadri-Momentum of the ejectile in the cm
  TLorentzVector ResidualMomentumEnergy(-EjectileMomentumCm,sqrt(pow(MomentumCmModule,2)+pow(Masses[1],2))); //Quadri-Momentum of the residual in the cm
  //Conversion of quantities in the lab frame
  EjectileMomentumEnergy.Boost(TotalMomentum.BoostVector()); //Boosting the ejectile cm momentum in the lab frame
  ResidualMomentumEnergy.Boost(TotalMomentum.BoostVector()); //Boosting the residual cm momentum in the lab frame
  //
  //Setting momenta of each particle
  fTheReactionProducts[0]->fMomentum=EjectileMomentumEnergy;
  fTheReactionProducts[1]->fMomentum=ResidualMomentumEnergy;
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
    } else SecondaryDecay(fTheReactionProducts[i],AnEvent,fTheAngularDistributions[i]->fSecondaryDecayAngularDistribution[0]); //If the particle has a decay pattern, process the decay pattern recursively till the end of the chain
  }
  //
  
  return AnEvent;
}

//The method SecondaryDecay takes care of calculating the secondary decay of 
//all particles originating by an initial particle in the exit channel and
//write all of them in the event used as the output of the GetEvent method.
//____________________________________________________
void UNISSequentialDecayTwoBody::SecondaryDecay(UNISIon * TheParticle, std::vector<UNISIon> & TheEvent, UNISAngularDistribution * ang_distr)
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
  //Generation of the secondary decay
  //Simulation of the primary reaction (2-body)
  TH1D * PrimaryDecayAngDistr = ang_distr ? ang_distr->fAngularDistribution : 0;  
  double theta_cms = PrimaryDecayAngDistr ? PrimaryDecayAngDistr->GetRandom()*TMath::DegToRad() : acos(gRandom->Uniform(0,1));
  double phi = gRandom->Uniform(0,2*TMath::Pi());
  const double InvariantMass=TheParticle->fMomentum.M();
  const double MomentumCmModule=sqrt((pow(InvariantMass,4)+pow((pow(Masses[0],2)-pow(Masses[1],2)),2)-2*pow(InvariantMass,2)*(pow(Masses[0],2)+pow(Masses[1],2)))/(4*pow(InvariantMass,2))); //Module of the momentum in the cm
  TVector3 EjectileMomentumCm(MomentumCmModule*sin(theta_cms)*cos(phi),MomentumCmModule*sin(theta_cms)*sin(phi),MomentumCmModule*cos(theta_cms)); //Momentum of the ejectile in the cm
  TLorentzVector EjectileMomentumEnergy(EjectileMomentumCm,sqrt(pow(MomentumCmModule,2)+pow(Masses[0],2))); //Quadri-Momentum of the ejectile in the cm
  TLorentzVector ResidualMomentumEnergy(-EjectileMomentumCm,sqrt(pow(MomentumCmModule,2)+pow(Masses[1],2))); //Quadri-Momentum of the residual in the cm
  //Conversion of quantities in the lab frame
  EjectileMomentumEnergy.Boost(TheParticle->fMomentum.BoostVector()); //Boosting the ejectile cm momentum in the lab frame
  ResidualMomentumEnergy.Boost(TheParticle->fMomentum.BoostVector()); //Boosting the residual cm momentum in the lab frame
  //
  //Setting momenta of each particle
  (TheParticle->fSecondaryParticles)[0]->fMomentum=EjectileMomentumEnergy;
  (TheParticle->fSecondaryParticles)[1]->fMomentum=ResidualMomentumEnergy;
  //
  
  //
  //Writing particles to the event or process additional secondary decays (recursively)
  for(int i=0; i<NumSecondaryParticles; i++)
  {
    if((TheParticle->fSecondaryParticles)[i]->fSecondaryParticles.size()==0) {
      TheEvent.push_back(*(TheParticle->fSecondaryParticles)[i]); //The particle does not decay, writing particle in the event
    } else SecondaryDecay((TheParticle->fSecondaryParticles)[i],TheEvent,fTheAngularDistributions[i]->fSecondaryDecayAngularDistribution[0]); //If the particle has a decay pattern, process the decay pattern recursively till the end of the chain
  }
  //
  
  return;
}

//____________________________________________________
int UNISSequentialDecayTwoBody::ProcessSetCommand(const char * line)
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
  } else if(WhatToSet.compare("ang_distr")==0) {

    //Setting the angular distribution of a given decay
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
    std::string AngularDistributionFile;
    LineStream>>AngularDistributionFile;
    AngularDistributionFile.assign(AngularDistributionFile.substr(AngularDistributionFile.find("\"")+1,AngularDistributionFile.find_last_of("\"")-(AngularDistributionFile.find("\"")+1)));
    //
    //Adding the new angular distribution
    TH1D * NewAngDistr = new TH1D(Form("AngDistr_%s",ANucleus.c_str()),"",500,0,180);
    const int lines_read=ReadAngularDistribution(AngularDistributionFile.c_str(),NewAngDistr);
    //
    //Setting quantities to the correct nucleus
    UNISAngularDistribution * TheAngDistrToSet=fTheAngularDistributions[PathToTheNucleus[0]];
    for(int i=1; i<PathToTheNucleus.size(); i++) {
      if(TheAngDistrToSet->fSecondaryDecayAngularDistribution.find(PathToTheNucleus[i])==TheAngDistrToSet->fSecondaryDecayAngularDistribution.end()) {
        TheAngDistrToSet->fSecondaryDecayAngularDistribution[PathToTheNucleus[i]]=new UNISAngularDistribution;
        TheAngDistrToSet->fSecondaryDecayAngularDistribution[PathToTheNucleus[i]]->fAngularDistribution=0;
      }
      TheAngDistrToSet=TheAngDistrToSet->fSecondaryDecayAngularDistribution[PathToTheNucleus[i]];
    }
    //
    TheAngDistrToSet->fAngularDistribution=lines_read>0 ? NewAngDistr : 0;
  } else {
    return 0; 
  }
  
  return 1;
}

//____________________________________________________
int UNISSequentialDecayTwoBody::ProcessDefineCommand(const char * line)
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
    fTheAngularDistributions[ParticleToSet]=new UNISAngularDistribution();
    fTheAngularDistributions[ParticleToSet]->fAngularDistribution=0;
    //
  } else {
    return 0; 
  }

  return 1;
}

//____________________________________________________
int UNISSequentialDecayTwoBody::ReadAngularDistribution(const char * file_name, TH1D * TheHisto)
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
    
    double x;
    double y;
    LineStream>>x>>y;
    
    TheHisto->Fill(x,y);
    
    NRead++;
  }
  
  return NRead;
}
