#ifndef TDETECTIONSETUP_H
#define TDETECTIONSETUP_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <TROOT.h>
#include <TEveManager.h>
#include <TGLViewer.h>
#include <TGFrame.h>

#include "TDetectionUnit.h"

class TDetectionSetup
{
public :
  TDetectionSetup(const char *);
  ~TDetectionSetup();
  
  void RegisterUnit(TDetectionUnit *);
  int Size() const;
  const char * GetName() const;
  
  bool IsInside(double theta, double phi, double x0=0, double y0=0, double z0=0) const;
  int GetDetectorIndex(double theta, double phi, double x0=0, double y0=0, double z0=0) const;
  TDetectionUnit * GetDetector(double theta, double phi, double x0=0, double y0=0, double z0=0); // returns the pointer to the detector object that detects the particle, 0 if the particle is not inside the cluster
  TDetectionUnit * GetDetector(int);  // Get a pointer to the detector at a certain index in the setup
  
  void Draw (Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const; // Draw the whole cluster on the x-y plane
  void Draw3D(Option_t * opt="") const; // Draw the whole cluster on a 3D space
  
private :
  std::string fName;
  int fNumDetectors;
  std::vector <TDetectionUnit*> fTheDetectors;
  
};

#endif
