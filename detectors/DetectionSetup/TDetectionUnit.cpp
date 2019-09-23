#include "TDetectionUnit.h"

//____________________________________________________
int TDetectionUnit::IsInside(double theta, double phi, double x0, double y0, double z0)
{
  return 0;
}

//____________________________________________________
int TDetectionUnit::GetPixel(double theta, double phi, double x0, double y0, double z0)
{
  return -1;
}

//____________________________________________________
TVector3 TDetectionUnit::GetImpactPointLab(double theta, double phi, double x0, double y0, double z0)
{
  return TVector3(0,0,0);
}

//____________________________________________________
void TDetectionUnit::Draw (Option_t * opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return;
}

//____________________________________________________
void TDetectionUnit::Draw3D (Option_t * opt) const
{
  return;
}
