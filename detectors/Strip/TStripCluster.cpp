#include "TStripCluster.h"

//____________________________________________________
TStripCluster::TStripCluster() :
fNumDetectors(0)
{}

//____________________________________________________
TStripCluster::~TStripCluster()
{
  for(int i=0; i<fTheDetectors.size(); i++) {
    delete fTheDetectors[i];   
  }
  fTheDetectors.clear();
}

//____________________________________________________
void TStripCluster::Clear()
{
  for(int i=0; i<fTheDetectors.size(); i++) {
    delete fTheDetectors[i];   
  }
  fTheDetectors.clear();
  fNumDetectors=0;    
}

//____________________________________________________
void TStripCluster::AddDetector(Double_t distance, Double_t theta_pos, Double_t phi_pos,
                                Int_t N_Strips, Double_t strip_width, Double_t inter_width, 
                                Double_t frame_width, Double_t dead_layer, Option_t * opt)
{
  TStripDetector * NewDetector = new TStripDetector(distance, theta_pos, phi_pos, N_Strips, strip_width, inter_width, frame_width, dead_layer, opt);
  fTheDetectors.push_back(NewDetector);
  fNumDetectors++;
}

//____________________________________________________
void TStripCluster::AddDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X,
                                Int_t N_Strips, Double_t strip_width, Double_t inter_width, 
                                Double_t frame_width, Double_t dead_layer, Option_t * opt)
{
  TStripDetector * NewDetector = new TStripDetector(X0, Y0, Z0, tilt_X, N_Strips, strip_width, inter_width, frame_width, dead_layer, opt);
  fTheDetectors.push_back(NewDetector);
  fNumDetectors++;
}

//____________________________________________________
bool TStripCluster::CheckOverlap() const
{
  return false;
}

//____________________________________________________
int TStripCluster::Size() const
{
  return fNumDetectors;
}

//____________________________________________________
TStripDetector * TStripCluster::GetDetector(int numdetector)
{
  return fTheDetectors[numdetector];
}

//____________________________________________________
int TStripCluster::IsInside(double theta, double phi, double x0, double y0, double z0) const
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return 1;   
  }
  return 0;
}

//____________________________________________________
TStripDetector * TStripCluster::GetDetector(double theta, double phi, double x0, double y0, double z0)
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return fTheDetectors[i];   
  }
  return 0;
}

//____________________________________________________
int TStripCluster::GetDetectorIndex(double theta, double phi, double x0, double y0, double z0) const
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return i;   
  }
  return -1;
}

//____________________________________________________
int TStripCluster::GetDetectorStripFront(int numdet, double theta, double phi, double x0, double y0, double z0) const
{
  return fTheDetectors[numdet]->GetStripFront(theta, phi, x0, y0, z0);   
}

//____________________________________________________
int TStripCluster::GetDetectorStripBack(int numdet, double theta, double phi, double x0, double y0, double z0) const
{
  return fTheDetectors[numdet]->GetStripBack(theta, phi, x0, y0, z0);   
}

//____________________________________________________
double TStripCluster::GetDetectorPixelTheta(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetThetaPixel(stripf, stripb);
}

//____________________________________________________
double TStripCluster::GetDetectorPixelPhi(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetPhiPixel(stripf, stripb); 
}

//____________________________________________________
TVector3 TStripCluster::GetDetectorPixelCenter(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetPixelCenter(stripf, stripb); 
}

//____________________________________________________
TVector3 TStripCluster::GetDetectorCenter(int numdet) const
{
  return fTheDetectors[numdet]->GetDetectorCenter();
}

//____________________________________________________
void TStripCluster::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw("", Xmin, Xmax, Ymin, Ymax) : fTheDetectors[i]->Draw("SAME");
  }
  
  return; 
}

//____________________________________________________
void TStripCluster::Draw3D(Option_t * draw_opt) const
{
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw3D("") : fTheDetectors[i]->Draw3D("SAME");
  }
  
  return; 
}

#ifdef GRAPHICAL_DEBUG
//____________________________________________________
void TStripCluster::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0, double Xmin, double Xmax, double Ymin, double Ymax)
{
  Draw("", Xmin, Xmax, Ymin, Ymax); //Crate the canvas and draw the cluster
  
  if(IsInside(theta_inc,phi_inc, x0, y0, z0))
  { 
    TStripDetector * TheDetector = GetDetector(theta_inc,phi_inc,x0,y0,z0);
    
    if(TheDetector==0) return;
    
    int nstripf = TheDetector->GetStripFront(theta_inc,phi_inc,x0,y0,z0);
    int nstripb = TheDetector->GetStripBack(theta_inc,phi_inc,x0,y0,z0);
    TVector3 TheImpactPoint = TheDetector->GetImpactPointLab(theta_inc,phi_inc,x0,y0,z0);
					   
    TGraph * impact_point = new TGraph();
    impact_point->SetPoint(0,TheImpactPoint.Y(),TheImpactPoint.X());
    impact_point->SetMarkerColor(kRed);
    
    impact_point->Draw("* SAME");
    
    if(nstripf>=0 && nstripb>=0) {
      TVector3 TheReconstructedImpactPoint = TheDetector->GetPixelCenter(nstripf, nstripb);
      TGraph * reconstructed_impact_point = new TGraph();
      reconstructed_impact_point->SetPoint(0,TheReconstructedImpactPoint.Y(),TheReconstructedImpactPoint.X());
      reconstructed_impact_point->SetMarkerColor(kBlue);
      reconstructed_impact_point->Draw("* SAME");
    }
  }
}
#endif

