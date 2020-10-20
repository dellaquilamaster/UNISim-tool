/* *****************************************
 * 17/10/2020
 * Class to handle Single-Sided Silicon Strip Detectors 
 * like those of Micron Semiconductors.
 * The class is derived by the former UNISStripDetector
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISSTRIPSINGLESIDEDDETECTOR_H
#define UNISSTRIPSINGLESIDEDDETECTOR_H

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

#define GRAPHICAL_DEBUG

/* *****************************************
 * 
 * Example
 * 
 * 
 * The reference frame used for the calculations has the z as the beam axis. In the upstream view, 
 * x is vertical and y goes towards the right hand side.
 * By defayult, detectors at phi=0 are located on the horizontal plane at the right side of the beam (upstream) <- this is no longer true
 * *****************************************/

class UNISStripSingleSidedDetector : public UNISDetectionUnit
{
private: 
  TVector3   TXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   TYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   TZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3   TTelescopeCenter; /*Telescope's center TVector3*/
  TVector3   TXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   TYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   TTrueImpactPoint; /*True impact point coordinates in the laboratory frame*/
  TVector3   TTelescopeImpactPoint; /*Impact point coordinates in the teloescope frame*/
  TVector3   TTopLeftCorner;  /*Telescope Top Left Corner*/
  TVector3   TTopRightCorner; /*Telescope Top Right Corner*/
  TVector3   TBottomLeftCorner; /*Telescope Bottom Left Corner*/
  Double_t   TTiltXAngle; /*tilt angle with respect to the horizontal (X) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Double_t   TTiltYAngle; /*tilt angle with respect to the vertical (Y) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Int_t      TStrips_number; /*Total number of strips*/
  Double_t   TStripTrue_width; /*width of each pixel in mm*/ 
  Double_t   TStripTrue_semi; /*half-width of the pixel including the frame*/
  Double_t   TInter_width; /*Width of the inter-strip of each strip*/
  Double_t   TFrame_width; /*Width of the frame of the telescope*/
  Double_t   TDeadLayer; /*external layer of the silicon that is a dead region*/
  Double_t   TStripEffective_width; /*effective width of each pixel in mm*/ 
  Double_t   TStripEffective_semi; /*half-width of the effective area of each pixel in mm*/ 
  Double_t   TTelescopeEffective_semi; /*half widht of the telescope's effective area in mm*/
  Double_t   TTelescopeTrue_semi; /*half widht of the entire telescope in mm including the frame*/
  Double_t   TTopLeftXCorner; /*X coordinate of the Telescope's top left corner*/
  Double_t   TTopLeftYCorner; /*Y coordinate of the Telescope's top left corner*/
  Double_t   TImpactX; /*X coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   TImpactY; /*Y coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   TImpactXprime; /*X coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Double_t   TImpactYprime; /*Y coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Double_t   Ta; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   Tb;
  Double_t   Tc;
  Double_t   Td;
  
public:
  UNISStripSingleSidedDetector(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, 
                    Int_t N_Strips=16, Double_t strip_width=0.312, Double_t inter_width=0.01, Double_t frame_width=0.3, Double_t dead_layer=0., Option_t * opt="");
  UNISStripSingleSidedDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X=0., Double_t tilt_Y=0.,
                    Int_t N_Strips=16, Double_t strip_width=0.312, Double_t inter_width=0.01, Double_t frame_width=0.3, Double_t dead_layer=0., Option_t * opt="");     
  ~UNISStripSingleSidedDetector();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area
  void      RotateX(Double_t);
  void      RotateY(Double_t);
  void      RotateZ(Double_t);
  void      Translate(Double_t, Double_t, Double_t);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TVector3  GetDetectorCenter(); // returns a TVector3 representing the center of the detector in the lab reference frame
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
private :
  TGeoVolume * fFrame;
  TGeoVolume * fStripGround;
  TGeoVolume * fStrip;
  //
  void Generate3D();
  void Rotate3DX(Double_t);
  void Rotate3DY(Double_t);
  void Rotate3DZ(Double_t);
  void Translate3D(Double_t, Double_t, Double_t);
  
#ifdef GRAPHICAL_DEBUG
  void      ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
#endif
} ;

#endif
