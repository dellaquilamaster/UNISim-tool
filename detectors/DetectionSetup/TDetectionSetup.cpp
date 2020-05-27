#include "TDetectionSetup.h"

//____________________________________________________
TDetectionSetup::TDetectionSetup(const char * name) :
fName(name),
fNumDetectors(0)
{}

//____________________________________________________
TDetectionSetup::~TDetectionSetup()
{}

//____________________________________________________
void TDetectionSetup::RegisterUnit(TDetectionUnit * NewUnit)
{
  fTheDetectors.push_back(NewUnit);
  fNumDetectors++;
}

//____________________________________________________
int TDetectionSetup::Size() const
{
  return fNumDetectors;
}

//____________________________________________________
const char * TDetectionSetup::GetName() const
{
  return fName.c_str();
}

//____________________________________________________
void TDetectionSetup::Draw (Option_t * opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw(opt, Xmin, Xmax, Ymin, Ymax) : fTheDetectors[i]->Draw(Form("%s SAME",opt));
  }
  
  return;
}

//____________________________________________________
void TDetectionSetup::Draw3D (Option_t * opt) const
{  
  for(Int_t i=0; i<fNumDetectors; i++)
  {
    i==0 ? fTheDetectors[i]->Draw3D(opt) : fTheDetectors[i]->Draw3D(Form("%s SAME",opt));
  }
  
  //
  if(gEve) {
    TGLViewer *v = gEve->GetDefaultGLViewer();
    v->ColorSet().Background().SetColor(kWhite);
    v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
    v->RefreshPadEditor(v);
  }
  //
  
  return; 
}

//____________________________________________________
TDetectionUnit * TDetectionSetup::GetDetector(int numdetector)
{
  return (TDetectionUnit*)fTheDetectors[numdetector];
}

//____________________________________________________
bool TDetectionSetup::IsInside(double theta, double phi, double x0, double y0, double z0) const
{
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) return true;   
  }
  return false;
}

//____________________________________________________
TDetectionUnit * TDetectionSetup::GetDetector(double theta, double phi, double x0, double y0, double z0)
{
  //
  //WARNING: more than 1 detector might see a particle
  //To describe shadowing in a correct way, if more than 1 detectors
  //detect the same particle, the one closer to the beam-line is considered
  std::map<double, int> HitPositionDistance;
  //
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) {
      HitPositionDistance[(fTheDetectors[i]->GetImpactPointLab(theta, phi, x0, y0, z0)-TVector3(x0,y0,z0)).Mag()]=i;   
    }
  }
  //
  if(HitPositionDistance.empty()) return 0;
  //
  auto min_distance_hit = std::min_element(HitPositionDistance.begin(), HitPositionDistance.end());
  //
  
  return fTheDetectors[min_distance_hit->second];
}

//____________________________________________________
int TDetectionSetup::GetDetectorIndex(double theta, double phi, double x0, double y0, double z0) const
{
  //
  //WARNING: more than 1 detector might see a particle
  //To describe shadowing in a correct way, if more than 1 detectors
  //detect the same particle, the one closer to the beam-line is considered
  std::map<double, int> HitPositionDistance;
  //
  for(int i=0; i<fNumDetectors; i++) {
    if(fTheDetectors[i]->IsInside(theta, phi, x0, y0, z0)) {
      HitPositionDistance[(fTheDetectors[i]->GetImpactPointLab(theta, phi, x0, y0, z0)-TVector3(x0,y0,z0)).Mag()]=i;   
    }
  }
  //
  if(HitPositionDistance.empty()) return -1;
  //
  auto min_distance_hit = std::min_element(HitPositionDistance.begin(), HitPositionDistance.end());
  //
  
  return min_distance_hit->second;
}

