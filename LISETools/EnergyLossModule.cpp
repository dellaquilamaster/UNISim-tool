#include "EnergyLossModule.h"


//______________________________________________
EnergyLossModule::EnergyLossModule(const char * program_path) :
fProgramPath(program_path),
fEnergyLossPrecision(0.001)
{
  if(fProgramPath.empty()) {
    fProgramPath.assign(".");
  }
  if(fProgramPath.find_last_of('/')==fProgramPath.length()-1) {
    fProgramPath.assign(fProgramPath.substr(0,fProgramPath.find_last_of('/')));
  }

  NucData=new nuclear_masses(Form("%s/input/masses.conf", fProgramPath.c_str()));
}

//______________________________________________
EnergyLossModule::~EnergyLossModule()
{
  delete NucData;
}

//______________________________________________
void EnergyLossModule::Clear()
{
  for(std::map<std::string,std::vector <double>>::iterator it = ParticleEnergy.begin(); it!=ParticleEnergy.end(); it++) {
    std::string TheElement(it->first);
    (it->second).clear();
    for(int i=0; i<NUM_MODELS; i++) {
      if((LiseELoss[i])[TheElement].size()) {
        (LiseELoss[i])[TheElement].clear();
      }
    }
    ParticleEnergyRange[TheElement].clear();
    for(int i=0; i<NUM_RANGE_MODELS; i++) {
      if((LiseRange[i])[TheElement].size()) {
        (LiseRange[i])[TheElement].clear();
      }
    }
  }
}

//______________________________________________
int EnergyLossModule::LoadEnergyLossFile(const char * file_name)
{
  std::ifstream FileIn(file_name);
  if(!FileIn.is_open()) {
    return -1;
  }

  std::string TheElement (file_name);
  TheElement.assign(TheElement.substr(TheElement.find("LISE_ELoss"),TheElement.find_last_of(".dat")-3-TheElement.find("LISE_ELoss")));

  int NRead=0;

  while (!FileIn.eof())
  {
    std::string LineRead;
    std::getline(FileIn, LineRead);

    LineRead.assign(LineRead.substr(0,LineRead.find('*')));
    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineRead);

    double energy;
    double eloss;

    for(int i=0; i<NUM_MODELS; i++) {
      LineStream >> energy >> eloss;
      (LiseELoss[i])[TheElement].push_back(eloss);
    }

    ParticleEnergy[TheElement].push_back(energy);

    NRead++;
  }

  Emin[TheElement]=(ParticleEnergy[TheElement])[0];
  Emax[TheElement]=(ParticleEnergy[TheElement])[ParticleEnergy[TheElement].size()-1];

  for(int i=0; i<NUM_MODELS; i++) {
    (SplineInterpolator[i])[TheElement] = new ROOT::Math::Interpolator();
    (SplineInterpolator[i])[TheElement]->SetData(ParticleEnergy[TheElement],LiseELoss[i][TheElement]);
  }

  return NRead;
}

//______________________________________________
int EnergyLossModule::LoadRangeFile(const char * file_name)
{
  std::ifstream FileIn(file_name);
  if(!FileIn.is_open()) {
    return -1;
  }

  std::string TheElement (file_name);
  TheElement.assign(TheElement.substr(TheElement.find("LISE_Range"),TheElement.find_last_of(".dat")-3-TheElement.find("LISE_Range")));

  int NRead=0;

  while (!FileIn.eof())
  {
    std::string LineRead;
    std::getline(FileIn, LineRead);

    if(LineRead.empty()) continue;
    LineRead.assign(LineRead.substr(0,LineRead.find('*')));
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineRead);

    double energy;
    double range;

    for(int i=0; i<NUM_RANGE_MODELS; i++) {
      LineStream >> energy >> range;
      (LiseRange[i])[TheElement].push_back(range);
    }

    ParticleEnergyRange[TheElement].push_back(energy);

    NRead++;
  }

  ERangeMin[TheElement]=(ParticleEnergyRange[TheElement])[0];
  ERangeMax[TheElement]=(ParticleEnergyRange[TheElement])[ParticleEnergyRange[TheElement].size()-1];

  for(int i=0; i<NUM_RANGE_MODELS; i++) {
    (RangeSplineInterpolator[i])[TheElement] = new ROOT::Math::Interpolator();
    (EnergyFromRangeSplineInterpolator[i])[TheElement] = new ROOT::Math::Interpolator();
    (RangeSplineInterpolator[i])[TheElement]->SetData(ParticleEnergyRange[TheElement],(LiseRange[i])[TheElement]);
    (EnergyFromRangeSplineInterpolator[i])[TheElement]->SetData((LiseRange[i])[TheElement],ParticleEnergyRange[TheElement]);
    (RangeMin[i])[TheElement]=(LiseRange[i])[TheElement][0];
    (RangeMax[i])[TheElement]=(LiseRange[i])[TheElement][(LiseRange[i])[TheElement].size()-1];
  }

  return NRead;
}

//______________________________________________
void EnergyLossModule::SetEnergyLossPrecision(double new_precision)
{
  fEnergyLossPrecision=new_precision;
} 

//______________________________________________
double EnergyLossModule::GetEnergyLoss(int Z, int A, double Einc, const char * material, double thickness_um, int model)
{
  double dThicknessMax=thickness_um; //This is the maximum thickness used as the integration step
  double IntegrateThickness=0;
  double dThickness=dThicknessMax;
  double Eresidual=Einc;
  double ELoss=0;

  double mass_uma=NucData->get_mass_Z_A_uma(Z,A);

  std::string TheCombination (Form("LISE_ELoss_Z%02d_A%02d_%s",Z, A, material));

  //Check if information is present
  if(SplineInterpolator[model].find(TheCombination)==SplineInterpolator[model].end()) {
    if(LoadEnergyLossFile(Form("%s/input/LISE_ELoss_Z%02d_A%02d_%s.dat", fProgramPath.c_str(), Z, A, material))<=0) {
      printf("Error: information not present for Z=%d A=%d material=%s\n", Z, A, material);
      return -100;
    }
 }

 //Retrieve the LISE information
 ROOT::Math::Interpolator * TheSplineInterpolator = (SplineInterpolator[model])[TheCombination];
 double TheEmin = Emin[TheCombination];

  for(;IntegrateThickness<thickness_um; IntegrateThickness+=dThickness)
  {
    if(Eresidual<=TheEmin*mass_uma) { //the particle stopped in the material
      ELoss=Einc;
      return ELoss;
    }

    if(TheSplineInterpolator->Deriv(Eresidual/mass_uma)!=0) {
      dThickness=fmin(dThicknessMax,std::fabs(fEnergyLossPrecision/TheSplineInterpolator->Deriv(Eresidual/mass_uma))); //variable integration step with fixed precision
    }

    double ELossStep=dThickness*TheSplineInterpolator->Eval(Eresidual/mass_uma);

    ELoss+=ELossStep;
    Eresidual-=ELossStep;
  }

  return ELoss;
}

//______________________________________________
double EnergyLossModule::GetResidualEnergy(int Z, int A, double Eloss, const char * material, double thickness_um, int model)
{
  double Einc=GetIncidentEnergy(Z,A,Eloss,material,thickness_um, model);
  if(Einc<0) return -1; //the particle cannot deposit this energy (energy greater than punch through energy)
  return Einc-Eloss;
}

//______________________________________________
double EnergyLossModule::GetIncidentEnergy(int Z, int A, double Eloss, const char * material, double thickness_um, int model)
{
  double EincStep=Eloss;
  double ElossStep;
  double dE=30.;

  ElossStep=GetEnergyLoss(Z,A,EincStep,material,thickness_um, model);

  if(ElossStep<Eloss) return -1; //the particle cannot deposit this energy (energy greater than punch through energy)

  for(;;EincStep+=dE)
  {
    ElossStep=GetEnergyLoss(Z,A,EincStep,material,thickness_um, model);

    if(ElossStep<Eloss) {
      dE=-std::fabs(dE)/2;
    }
    if(ElossStep>Eloss && dE<0) {
      dE=std::fabs(dE);
    }
    if(std::fabs(dE)<fEnergyLossPrecision) break;
  }

  return EincStep;
}

//______________________________________________
double EnergyLossModule::GetIncidentEnergyFromResidual(int Z, int A, double Eres, const char * material, double thickness_um, int model)
{
  double EincStep=Eres;
  double E_residual_step;
  double dE=30.;

  while(1)
  {
    double EincTest = EincStep+dE;
    E_residual_step=EincTest-GetEnergyLoss(Z,A,EincTest,material,thickness_um, model);
    
    if(E_residual_step<Eres) {
      EincStep=EincTest;
      dE=1.1*dE;
    } else {
      dE=0.1*dE; 
    }
    if(dE<fEnergyLossPrecision) break;
  }

  return EincStep;
}

//______________________________________________
double EnergyLossModule::GetThicknessFromEnergyLoss(int Z, int A, double Einc, double Eloss, const char * material, int model)
{
  //
  if(Einc<Eloss) return -1; //The particle does not have sufficient energy to loose Eloss
  //
  double Range=GetRangeFromEnergy(Z,A,Einc,material,model);
  if(Range<0) return -100;
  //
  double thickness_precision_target = 0.01;
  double MinThickness=0.;
  double MaxThickness=Range;
  double TestThickness=(MinThickness+MaxThickness)/2.;
  double ElossStep=Einc;
  //
  
  while(1)
  {
    //
    if(ElossStep<Eloss) {
      //Test thickness is smaller than the target value
      MinThickness=TestThickness;
    } else {
      //Test thickness is greater than the target value
      MaxThickness=TestThickness;
    }
    //
    
    //
    TestThickness=(MinThickness+MaxThickness)/2.;
    //
    
    //
    if(std::fabs(MaxThickness-MinThickness)<thickness_precision_target) break;
    //
    
    ElossStep=GetEnergyLoss(Z,A,Einc,material,TestThickness, model);

  }

  return TestThickness;
}

//______________________________________________
double EnergyLossModule::GetRangeFromEnergy(int Z, int A, double Einc, const char * material, int model)
{
  double mass_uma=NucData->get_mass_Z_A_uma(Z,A);
  std::string TheCombination (Form("LISE_Range_Z%02d_A%02d_%s",Z, A, material));

  //Check if information is present
  if(RangeSplineInterpolator[model].find(TheCombination)==RangeSplineInterpolator[model].end()) {
    if(LoadRangeFile(Form("%s/input/LISE_Range_Z%02d_A%02d_%s.dat", fProgramPath.c_str(), Z, A, material))<=0) {
      printf("Error: range information not present for Z=%d A=%d material=%s\n", Z, A, material);
      return -100;
    }
  }

  //Retrieve the LISE information
  ROOT::Math::Interpolator * TheRangeSplineInterpolator = (RangeSplineInterpolator[model])[TheCombination];
  double TheERangeMin = ERangeMin[TheCombination];
  double TheERangeMax = ERangeMax[TheCombination];

  if(Einc/mass_uma>=TheERangeMin && Einc/mass_uma<=TheERangeMax) {
    return TheRangeSplineInterpolator->Eval(Einc/mass_uma);
  } else return -1;
}

//______________________________________________
double EnergyLossModule::GetEnergyFromRange(int Z, int A, double range, const char * material, int model)
{
  double mass_uma=NucData->get_mass_Z_A_uma(Z,A);
  std::string TheCombination (Form("LISE_Range_Z%02d_A%02d_%s",Z, A, material));

  //Check if information is present
  if(RangeSplineInterpolator[model].find(TheCombination)==RangeSplineInterpolator[model].end()) {
    if(LoadRangeFile(Form("%s/input/LISE_Range_Z%02d_A%02d_%s.dat", fProgramPath.c_str(), Z, A, material))<=0) {
      printf("Error: range information not present for Z=%d A=%d material=%s\n", Z, A, material);
      return -100;
    }
  }

  //Retrieve the LISE information
  ROOT::Math::Interpolator * TheEnergyFromRangeSplineInterpolator = (EnergyFromRangeSplineInterpolator[model])[TheCombination];
  double TheRangeMin = (RangeMin[model])[TheCombination];
  double TheRangeMax = (RangeMax[model])[TheCombination];

  if(range>=TheRangeMin && range<=TheRangeMax) {
    return mass_uma*TheEnergyFromRangeSplineInterpolator->Eval(range);
  } else return -1;
}

//______________________________________________
void EnergyLossModule::DrawdEdx(int Z, int A, const char * material, int model)
{
  std::string TheCombination (Form("LISE_ELoss_Z%02d_A%02d_%s",Z, A, material));

  //Check if information is present
  if(SplineInterpolator[model].find(TheCombination)==SplineInterpolator[model].end()) {
    if(LoadEnergyLossFile(Form("%s/input/LISE_ELoss_Z%02d_A%02d_%s.dat", fProgramPath.c_str(), Z, A, material))<=0) {
      printf("Error: information not present for Z=%d A=%d material=%s\n", Z, A, material);
      return;
    }
 }

  const int NPoints = ParticleEnergy[TheCombination].size();
  double E_LISE_Values[NPoints];
  double LISE_Values[NPoints];

  const int NPointsInterpolation = NPoints*1000;
  double E_LISE_ValuesInterpolation[NPointsInterpolation];
  double LISE_ValuesInterpolation[NPointsInterpolation];

  for(int i=0; i<NPoints; i++) {
    E_LISE_Values[i]=(ParticleEnergy[TheCombination])[i];
    LISE_Values[i]=(LiseELoss[model])[TheCombination][i];
  }

  for(int i=0; i<NPointsInterpolation; i++) {
    E_LISE_ValuesInterpolation[i]=i*(Emax[TheCombination]-Emin[TheCombination])/NPointsInterpolation;
    LISE_ValuesInterpolation[i]=(SplineInterpolator[model])[TheCombination]->Eval(E_LISE_ValuesInterpolation[i]);
  }

  TGraph *LISEGraph = new TGraph(NPoints,E_LISE_Values,LISE_Values);
  TGraph *LISEGraphInterpolation = new TGraph(NPointsInterpolation,E_LISE_ValuesInterpolation,LISE_ValuesInterpolation);

  LISEGraph->Draw("A*");
  LISEGraphInterpolation->Draw("same L");
  LISEGraphInterpolation->SetLineColor(kRed);

  return;
}
