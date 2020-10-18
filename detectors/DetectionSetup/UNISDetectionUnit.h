#ifndef UNISDETECTIONUNIT_H
#define UNISDETECTIONUNIT_H

#include <TROOT.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>

class UNISDetectionUnit
{
public : 
  
  virtual int IsInside(double theta, double phi, double x0=0, double y0=0, double z0=0);
  virtual int GetPixel(double theta, double phi, double x0=0, double y0=0, double z0=0);
  virtual TVector3 GetImpactPointLab(double theta, double phi, double x0=0, double y0=0, double z0=0);

  virtual void Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const;
  virtual void Draw3D(Option_t * opt="") const;
  
protected :
  TGeoVolumeAssembly * fDetector;
  TGeoHMatrix * fDetectorMatrix;
  
};

#endif
