/* *****************************************
 * 18/10/2020
 * Class to handle OSCAR
 * The class utilizes the objects UNISSiliconPhotoDiode and UNISStripSingleSidedDetector
 * to compose an OSCAR telescope.
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISOSCARTELESCOPE_H
#define UNISOSCARTELESCOPE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <TVector3.h>
#include <TMath.h>
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
#include "../Strip/UNISStripSingleSidedDetector.h"
#include "../Silicon/UNISSiliconPhotoDiode.h"

/*
 * NOTE on the reference frame.
 * In the reference frame used to build OSCAR, z is the beam axis, x is vertical towards the top and y is horizontal towards the right.
 * Positioning the detector at phi=90 deg corresponds to the right side of the beamline looking downstream.
 * To look at the impact point downstream, X is towards left, and Y is towards the top
 * 
 */

class UNISOscarTelescope : public UNISDetectionUnit
{
private : 
  const int fNumPads;
  const int fNumStrips;
  const int fPadRowsColumns;
  const double fPadWidth;
  const double fCollimatorWidth;
  const double fPadSemi;
  const double fCollimatorSemi;
  const double fPadFrameWidth;
  const double fPadBottomContactsWidth;
  const double fPhotoDiodeWidth;
  const double fPhotoDiodeHeight;
  const double fFrameWidth;
  const double fFrameHeight;
  bool fIsStrip;
  bool fIsCollimator;
  const double fStripPadDistance;
  //
  TVector3   fXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   fYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   fZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3   fCenter; /*TVector3 representing the center of the active area*/
  TVector3   fXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   fYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   fLabImpactPoint;
  TVector3   fFrameImpactPoint;
  TVector3   fTopLeftCorner;  /*Telescope Top Left Corner*/
  TVector3   fTopRightCorner; /*Telescope Top Right Corner*/
  TVector3   fBottomLeftCorner; /*Telescope Bottom Left Corner*/
  TVector3   fBottomRightCorner; /*Telescope Bottom Right Corner*/
  Double_t   fa; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   fb;
  Double_t   fc;
  Double_t   fd;  
  //
  UNISSiliconPhotoDiode ** fPads;
  UNISStripSingleSidedDetector * fStrip;
  //
  
public :
  UNISOscarTelescope(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, Option_t * opt=""); //Creates an OSCAR telescope if the opt="", Creates only the second stage if opt="pads", Creates also pad collimators if the option "col" is specified
  UNISOscarTelescope(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y, Option_t * opt="");  //Creates an OSCAR telescope if the opt="", Creates only the second stage if opt="pads", Creates also pad collimators if the option "col" is specified   
  ~UNISOscarTelescope();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area
  Int_t     GetStrip(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
  Int_t     GetPad(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
  void      RotateX(Double_t);
  void      RotateY(Double_t);
  void      RotateZ(Double_t);
  void      Translate(Double_t, Double_t, Double_t);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TVector3  GetDetectorCenter(); // returns a TVector3 representing the center of the detector in the lab reference frame
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
private :
  TGeoVolume * fFramePads;
  TGeoVolume * fFramePreAmps;
  TGeoVolume * fPreAmps;
  TGeoVolume * fCollimatorFrame;
  TGeoVolume ** fCollimatorHole;
  TGeoVolume * fCollimator;
  //
  void Generate3D();
  
} ;

#endif
