/* *****************************************
 * 21/05/2020
 * Class to handle a simple Fazia quartet
 * facing the target perpendicularly at a distance 
 * of 100 cm.
 * The class is derived from TStrip
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISFAZIAQUARTET_H
#define UNISFAZIAQUARTET_H

#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TLine.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TEveManager.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoPara.h>
#include <TEveGeoShape.h>
#include <TEveTrans.h>
#include <TEveText.h>
#include <TGLViewer.h>
#include <TCanvas.h>

#include "../DetectionSetup/TDetectionUnit.h"

#define GRAPHICAL_DEBUG

/* *****************************************
 * 
 * Example
 * 
 * 2 x 2 square pads

 *  __________
 * | 0  | 1  | 
 * |----------
 * | 2  | 3  |   
 * |---------|
 * 
 * The reference frame used for the calculations has the z as the beam axis. In the upstream view, 
 * x is vertical and y goes towards the right hand side.
 * By defayult, detectors at phi=0 are located on the horizontal plane at the right side of the beam (upstream)
 * *****************************************/

class UNISFaziaQuartet : public TDetectionUnit
{
private: 
  TVector3    TXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3    TYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3    TZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3**  TCenters; /*Matrix containing pads's centers*/
  TVector3    TTelescopeCenter; /*Telescope's center TVector3*/
  TVector3    TXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3    TYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3    TTrueImpactPoint; /*True impact point coordinates in the laboratory frame*/
  TVector3    TTelescopeImpactPoint; /*Impact point coordinates in the teloescope frame*/
  TVector3    TTopLeftCorner;  /*Telescope Top Left Corner*/
  TVector3    TTopRightCorner; /*Telescope Top Right Corner*/
  TVector3    TBottomLeftCorner; /*Telescope Bottom Left Corner*/
  const Int_t TPads_number; /*Total number of pads*/
  const Int_t TRowColumn; /*total number of pads row or column*/
  Double_t    TNominalDistance; //distance between the center of the quartet and the target (each quartet is facing perpendicularly the target when placed at this nominal distance) - 100 cm
  Double_t    TDisplacement; //displacement from the nominal 100cm distance
  Double_t    TNominalTheta; //Theta angle identifying the center of the quartet when the array is placed at 100 cm from the target
  Double_t    TNominalPhi; //Phi angle identifying the center of the quartet when the array is placed at 100 cm from the target
  Double_t    TPadEffective_width; /*effective width of each pad in mm*/ 
  Double_t    TPadEffective_semi; /*half-width of the effective area of each pad in mm*/ 
  Double_t    TFrame_width; /*Width of the frame of the telescope*/
  Double_t    TPadTrue_width; /*width of each pad in mm*/ 
  Double_t    TPadTrue_semi; /*half-width of the pad including the frame*/
  Double_t    TTelescopeTrue_semi; /*half widht of the entire telescope in mm including the frame*/
  Double_t    TTopLeftXCorner; /*X coordinate of the Telescope's top left corner*/
  Double_t    TTopLeftYCorner; /*Y coordinate of the Telescope's top left corner*/
  Double_t**  TCentersXprime; /*Matrix containing pads's centers X corrdinates in the telescope frame respect to the top left corner*/
  Double_t**  TCentersYprime; /*Matrix containing pads's centers Y coordinates in the telescope frame respect to the top left corner*/  
  Double_t    TImpactX; /*X coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t    TImpactY; /*Y coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t    TImpactXprime; /*X coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Double_t    TImpactYprime; /*Y coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Double_t    Ta; /*parameters for the equation idetifying the plane of the detector*/
  Double_t    Tb;
  Double_t    Tc;
  Double_t    Td;
  
public:
  UNISFaziaQuartet(Double_t theta_pos=0, Double_t phi_pos=0, Double_t displacement=0., Double_t pad_width=2., Double_t frame_width=0.1);
  ~UNISFaziaQuartet();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pad fired, -1 if not inside the active area
  void      RotateX(Double_t);
  void      RotateY(Double_t);
  void      RotateZ(Double_t);
  void      Translate(Double_t, Double_t, Double_t);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TGraph*   GetGraphObject();
  TVector3  GetDetectorCenter(); // returns a TVector3 representing the center of the detector in the lab reference frame
  TVector3  GetPadCenter(int pad); // returns a TVector3 representing the center of the pad, identified by a given index, in the lab reference frame
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
#ifdef GRAPHICAL_DEBUG
  void      ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
#endif
  
private :
  TEveGeoShape * DetectorQuartet;
  TEveGeoShape *** fPad;
  TEveGeoShape *** fTopFrame;
  TEveGeoShape *** fBottomFrame;
  TEveGeoShape *** fLeftFrame;
  TEveGeoShape *** fRightFrame;
  TEveGeoShape *** fCsICrystal;
  //
  void Generate3D(Double_t theta, Double_t phi);
  void Rotate3DX(Double_t);
  void Rotate3DY(Double_t);
  void Rotate3DZ(Double_t);
  void Translate3D(Double_t, Double_t, Double_t);
} ;

#endif
