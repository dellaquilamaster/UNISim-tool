/* *****************************************
 * 18/10/2020
 * Class to handle a silicon PhotoDiode hamamatsu-like 
 * The class is derived by the former UNISStripSingleSidedDetector
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISCOLLIMATEDSILICON_H
#define UNISCOLLIMATEDSILICON_H

#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TLine.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoPara.h>
#include <TGeoShape.h>
#include <TGeoCompositeShape.h>

#include "../DetectionSetup/UNISDetectionUnit.h"

class UNISCollimatedSilicon : public UNISDetectionUnit
{
private: 
  TVector3   fXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   fYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   fZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3   fCenter; /*TVector3 representing the center of the active area*/
  TVector3   fXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   fYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   fLabImpactPoint; /*Iimpact point coordinates in the laboratory frame*/
  TVector3   fPadImpactPoint; /*Impact point coordinates in the pad frame*/
  TVector3   fTopLeftReference;  /*A reference a the top left of the center*/
  TVector3   fTopRightReference; /*A reference a the top right of the center*/
  TVector3   fBottomLeftReference; /*A reference at the bottom left of the center*/
  TVector3   fBottomRightReference; /*A reference at the bottom right of the center*/
  Double_t   fTiltXAngle; /*tilt angle with respect to the horizontal (X) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fTiltYAngle; /*tilt angle with respect to the vertical (Y) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fTiltZAngle; /*tilt angle with respect to the beam (Z) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fPadSemi; /*half-width of the pad*/
  Double_t   fCollimatorInnerRadius; /*Inner radius of the collimator*/
  Double_t   fCollimatorOuterRadius; /*Outer radius of the collimator*/
  Double_t   fBottomContactsWidth; /*Width of the space reserved to the electrical contacts (this adds up with the bottom frame)*/
  Double_t   fImpactX; /*X coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   fImpactY; /*Y coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   fa; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   fb;
  Double_t   fc;
  Double_t   fd;
  
public:
  UNISCollimatedSilicon(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, Double_t collimator_inner_radius=0.2, Double_t collimator_outer_radius=0.5);
  UNISCollimatedSilicon(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y, Double_t tilt_Z, Double_t collimator_inner_radius, Double_t collimator_outer_radius);     
  ~UNISCollimatedSilicon();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area
  void      RotateX(Double_t);
  void      RotateY(Double_t);
  void      RotateZ(Double_t);
  void      Translate(Double_t, Double_t, Double_t);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TVector3  GetDetectorCenter(); // returns a TVector3 representing the center of the detector active area in the lab reference frame
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
private :
  TGeoVolume * fCollimator;
  TGeoVolume * fActiveArea;
  //
  void Generate3D();
  void Rotate3DX(Double_t);
  void Rotate3DY(Double_t);
  void Rotate3DZ(Double_t);
  void Translate3D(Double_t, Double_t, Double_t);
  
} ;

#endif
