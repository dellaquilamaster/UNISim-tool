#include "LISETools/nuclear_masses.h"

void debugger()
{
  TChain * Data = new TChain("t");
  
  Data->Add("UNIS_10000events.root");
 
  Data->SetMakeClass(1);
  int fmulti;
  int fZ[10];
  int fA[10];
  double fKinEnergy[10];
  double fTheta[10];
  double fPhi[10];
  Data->SetBranchAddress("cluster.fmulti",&fmulti);
  Data->SetBranchAddress("cluster.fZ",fZ);
  Data->SetBranchAddress("cluster.fA",fA);
  Data->SetBranchAddress("cluster.fKinEnergy",fKinEnergy);
  Data->SetBranchAddress("cluster.fThetaOrigin",fTheta);
  Data->SetBranchAddress("cluster.fPhiOrigin",fPhi);
  Data->SetBranchStatus("*", false);
  Data->SetBranchStatus("cluster.fmulti", true);
  Data->SetBranchStatus("cluster.fZ", true);
  Data->SetBranchStatus("cluster.fA", true);
  Data->SetBranchStatus("cluster.fKinEnergy", true);
  Data->SetBranchStatus("cluster.fThetaOrigin", true);
  Data->SetBranchStatus("cluster.fPhiOrigin", true);
  
  Long64_t nentries = Data->GetEntries();
  
  nuclear_masses NucData("./LISETools/input/masses.conf");
  
  TH1D * HistoEx = new TH1D ("HistoEx","",300,0,30);
  
  for (Long64_t jentry = 0; jentry<nentries; jentry++)
  {
    Data->GetEntry(jentry);
    
    TLorentzVector Momento[10];
    
    for(int i =0; i<fmulti; i++) 
    {
      double massa = NucData.get_mass_Z_A(fZ[i],fA[i]);
      double P = sqrt(pow(fKinEnergy[i]+massa,2)-pow(massa,2));
      
      Momento[i]=TLorentzVector(P*sin(fTheta[i])*cos(fPhi[i]),P*sin(fTheta[i])*sin(fPhi[i]),P*cos(fTheta[i]),massa+fKinEnergy[i]);
    }
    
    for(int i=0; i<fmulti; i++) {
      if(fZ[i]==2 && fA[i]==6) {
        for(int j=i+1; j<fmulti; j++) {
          if(fZ[j]==2 && fA[j]==6) {
            if(fKinEnergy[i]>0 && fKinEnergy[j]>0)
              HistoEx->Fill((Momento[i]+Momento[j]).M()-NucData.get_mass_Z_A(fZ[i]+fZ[j],fA[i]+fA[j])); 
          }
        }
      }
    }
    
  }
  
  HistoEx->Draw();
  
  return;
}
