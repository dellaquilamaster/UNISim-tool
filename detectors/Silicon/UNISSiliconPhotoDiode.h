/* *****************************************
 * 18/10/2020
 * Class to handle a silicon PhotoDiode hamamatsu-like 
 * The class is derived by the former UNISStripSingleSidedDetector
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISSILICONPHOTODIODE_H
#define UNISSILICONPHOTODIODE_H

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

class UNISSiliconPhotoDiode : public UNISDetectionUnit
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
  TVector3   fTopLeftCorner;  /*Telescope Top Left Corner*/
  TVector3   fTopRightCorner; /*Telescope Top Right Corner*/
  TVector3   fBottomLeftCorner; /*Telescope Bottom Left Corner*/
  TVector3   fBottomRightCorner; /*Telescope Bottom Right Corner*/
  Double_t   fTiltXAngle; /*tilt angle with respect to the horizontal (X) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fTiltYAngle; /*tilt angle with respect to the vertical (Y) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fTiltZAngle; /*tilt angle with respect to the beam (Z) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   fPadWidth; /*width of the pad in cm*/ 
  Double_t   fPadEffectiveWidth; /*width of the active area of the pad in cm*/ 
  Double_t   fPadSemi; /*half-width of the pad*/
  Double_t   fPadEffectiveSemi; /*half-width of the effective area of the pad*/
  Double_t   fFrameWidth; /*Width of the frame*/
  Double_t   fBottomContactsWidth; /*Width of the space reserved to the electrical contacts (this adds up with the bottom frame)*/
  Double_t   fImpactX; /*X coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   fImpactY; /*Y coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   fa; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   fb;
  Double_t   fc;
  Double_t   fd;
  
public:
  UNISSiliconPhotoDiode(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, Double_t collimator_size=1.0);
  UNISSiliconPhotoDiode(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y, Double_t tilt_Z, Double_t collimator_size=1.0);     
  ~UNISSiliconPhotoDiode();

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
  TGeoVolume * fFrameLateral;
  TGeoVolume * fFrameHorizontal;
  TGeoVolume * fBottomContacts;
  TGeoVolume * fPad;
  //
  void Generate3D();
  void Rotate3DX(Double_t);
  void Rotate3DY(Double_t);
  void Rotate3DZ(Double_t);
  void Translate3D(Double_t, Double_t, Double_t);
  
} ;

#endif
