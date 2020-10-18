#include "UNISDetectionUnit.h"

//____________________________________________________
int UNISDetectionUnit::IsInside(double theta, double phi, double x0, double y0, double z0)
{
  return 0;
}

//____________________________________________________
int UNISDetectionUnit::GetPixel(double theta, double phi, double x0, double y0, double z0)
{
  return -1;
}

//____________________________________________________
TVector3 UNISDetectionUnit::GetImpactPointLab(double theta, double phi, double x0, double y0, double z0)
{
  return TVector3(0,0,0);
}

//____________________________________________________
void UNISDetectionUnit::Draw (Option_t * opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return;
}

//____________________________________________________
void UNISDetectionUnit::Draw3D (Option_t * opt) const
{
  return;
}
