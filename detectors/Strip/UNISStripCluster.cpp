#include "UNISStripCluster.h"

//____________________________________________________
UNISStripCluster::UNISStripCluster() :
fNumDetectors(0)
{}

//____________________________________________________
UNISStripCluster::~UNISStripCluster()
{
  for(unsigned int i=0; i<fTheDetectors.size(); i++) {
    delete fTheDetectors[i];   
  }
  fTheDetectors.clear();
}

//____________________________________________________
void UNISStripCluster::Clear()
{
  for(unsigned int i=0; i<fTheDetectors.size(); i++) {
    delete fTheDetectors[i];   
  }
  fTheDetectors.clear();
  fNumDetectors=0;    
}

//____________________________________________________
void UNISStripCluster::AddDetector(Double_t distance, Double_t theta_pos, Double_t phi_pos,
                                Int_t N_Strips, Double_t strip_width, Double_t inter_width, 
                                Double_t frame_width, Double_t dead_layer, Option_t * opt)
{
  UNISStripDetector * NewDetector = new UNISStripDetector(distance, theta_pos, phi_pos, N_Strips, strip_width, inter_width, frame_width, dead_layer, opt);
  fTheDetectors.push_back(NewDetector);
  fNumDetectors++;
}

//____________________________________________________
void UNISStripCluster::AddDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y,
                                Int_t N_Strips, Double_t strip_width, Double_t inter_width, 
                                Double_t frame_width, Double_t dead_layer, Option_t * opt)
{
  UNISStripDetector * NewDetector = new UNISStripDetector(X0, Y0, Z0, tilt_X, tilt_Y, N_Strips, strip_width, inter_width, frame_width, dead_layer, opt);
  fTheDetectors.push_back(NewDetector);
  fNumDetectors++;
}

//____________________________________________________
bool UNISStripCluster::CheckOverlap() const
{
  return false;
}

//____________________________________________________
int UNISStripCluster::Size() const
{
  return fNumDetectors;
}

//____________________________________________________
UNISStripDetector * UNISStripCluster::GetDetector(int numdetector)
{
  return fTheDetectors[numdetector];
}

//____________________________________________________
int UNISStripCluster::IsInside(double theta, double phi, double x0, double y0, double z0) const
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return 1;   
  }
  return 0;
}

//____________________________________________________
UNISStripDetector * UNISStripCluster::GetDetector(double theta, double phi, double x0, double y0, double z0)
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return fTheDetectors[i];   
  }
  return 0;
}

//____________________________________________________
int UNISStripCluster::GetDetectorIndex(double theta, double phi, double x0, double y0, double z0) const
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return i;   
  }
  return -1;
}

//____________________________________________________
int UNISStripCluster::GetDetectorStripFront(int numdet, double theta, double phi, double x0, double y0, double z0) const
{
  return fTheDetectors[numdet]->GetStripFront(theta, phi, x0, y0, z0);   
}

//____________________________________________________
int UNISStripCluster::GetDetectorStripBack(int numdet, double theta, double phi, double x0, double y0, double z0) const
{
  return fTheDetectors[numdet]->GetStripBack(theta, phi, x0, y0, z0);   
}

//____________________________________________________
double UNISStripCluster::GetDetectorPixelTheta(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetThetaPixel(stripf, stripb);
}

//____________________________________________________
double UNISStripCluster::GetDetectorPixelPhi(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetPhiPixel(stripf, stripb); 
}

//____________________________________________________
TVector3 UNISStripCluster::GetDetectorPixelCenter(int numdet, int stripf, int stripb) const
{
  return fTheDetectors[numdet]->GetPixelCenter(stripf, stripb); 
}

//____________________________________________________
TVector3 UNISStripCluster::GetDetectorCenter(int numdet) const
{
  return fTheDetectors[numdet]->GetDetectorCenter();
}

//____________________________________________________
void UNISStripCluster::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw("", Xmin, Xmax, Ymin, Ymax) : fTheDetectors[i]->Draw("SAME");
  }
  
  return; 
}

//____________________________________________________
void UNISStripCluster::Draw3D(Option_t * draw_opt) const
{
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw3D(draw_opt) : fTheDetectors[i]->Draw3D(Form("%s SAME", draw_opt));
  }
  
  return; 
}

#ifdef GRAPHICAL_DEBUG
//____________________________________________________
void UNISStripCluster::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0, double Xmin, double Xmax, double Ymin, double Ymax)
{
  Draw("", Xmin, Xmax, Ymin, Ymax); //Crate the canvas and draw the cluster
  
  if(IsInside(theta_inc,phi_inc, x0, y0, z0))
  { 
    UNISStripDetector * TheDetector = GetDetector(theta_inc,phi_inc,x0,y0,z0);
    
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

