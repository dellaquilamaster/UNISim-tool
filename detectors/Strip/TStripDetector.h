/* *****************************************
 * 22/12/2015
 * C++ Object used to create a simple 
 * Farcos-like strip detector made by a fixable 
 * number of pixels of fixable dimensions. 
 * Every cluster of pixels can be placed 
 * in any needed position.
 * 24/04/2019
 * Fully debugged. Geometry re-implemented.
 * New draw functions added.
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef _TStripDetector
#define _TStripDetector

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

#include "../DetectionSetup/TDetectionUnit.h"

#define GRAPHICAL_DEBUG

/* *****************************************
 * 
 * Example
 * 
 * 16 square pixels
 *   strips front are vertical
 *   0    1    2    3    <-- strips back are horizontal
 *  ___________________
 * | 0  | 1  | 2  | 3  |   0 back strip
 * |-------------------|
 * | 4  | 5  | 6  | 7  |   1 back strip
 * |-------------------|
 * | 8  | 9  | 10 | 11 |   2
 * |-------------------|
 * | 12 | 13 | 14 | 15 |   3
 * |___________________|
 * 
 * The reference frame used for the calculations has the z as the beam axis. In the upstream view, 
 * x is vertical and y goes towards the right hand side.
 * By defayult, detectors at phi=0 are located on the horizontal plane at the right side of the beam (upstream)
 * *****************************************/

class TStripDetector : public TDetectionUnit
{
private: 
  TVector3   TXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   TYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   TZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3** TCenters; /*Matrix containing pixels's centers*/
  TVector3   TTelescopeCenter; /*Telescope's center TVector3*/
  TVector3   TXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   TYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   TTrueImpactPoint; /*True impact point coordinates in the laboratory frame*/
  TVector3   TTelescopeImpactPoint; /*Impact point coordinates in the teloescope frame*/
  TVector3   TTopLeftCorner;  /*Telescope Top Left Corner*/
  TVector3   TTopRightCorner; /*Telescope Top Right Corner*/
  TVector3   TBottomLeftCorner; /*Telescope Bottom Left Corner*/
  Double_t   TTiltAngle; /*tilt angle with respect to the vertical (X) axis. This is -9999 is the detector is facing the target perpendicularly (first constructor)*/
  Int_t      TStrips_number; /*Total number of strips*/
  Int_t      TPixels_number; /*Total number of pixels*/
  Int_t      TRowColumn; /*total number of pixels row or column*/
  Double_t   TPixelTrue_width; /*width of each pixel in mm*/ 
  Double_t   TPixelTrue_semi; /*half-width of the pixel including the frame*/
  Double_t   TInter_width; /*Width of the inter-strip of each strip*/
  Double_t   TFrame_width; /*Width of the frame of the telescope*/
  Double_t   TDeadLayer; /*external layer of the silicon that is a dead region*/
  Double_t   TPixelEffective_width; /*effective width of each pixel in mm*/ 
  Double_t   TPixelEffective_semi; /*half-width of the effective area of each pixel in mm*/ 
  Double_t   TTelescopeEffective_semi; /*half widht of the telescope's effective area in mm*/
  Double_t   TTelescopeTrue_semi; /*half widht of the entire telescope in mm including the frame*/
  Double_t   TTopLeftXCorner; /*X coordinate of the Telescope's top left corner*/
  Double_t   TTopLeftYCorner; /*Y coordinate of the Telescope's top left corner*/
  Double_t** TCentersXprime; /*Matrix containing pixels's centers X corrdinates in the telescope frame respect to the top left corner*/
  Double_t** TCentersYprime; /*Matrix containing pixels's centers Y coordinates in the telescope frame respect to the top left corner*/  
  Double_t   TImpactX; /*X coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   TImpactY; /*Y coordinate of the impact point on the telescope surface in the telescope frame*/
  Double_t   TImpactXprime; /*X coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Double_t   TImpactYprime; /*Y coordinate of the impact point on the telescope surface in the telescope frame respect to the top left corner*/
  Int_t *    TStrip_hit; /*array containing strip front and back hit by the particle*/
  Double_t   Ta; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   Tb;
  Double_t   Tc;
  Double_t   Td;
  
public:
  TStripDetector(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, Int_t N_Strips=16, 
		   Double_t strip_width=0.2, Double_t inter_width=0.01, Double_t frame_width=0.2, Double_t dead_layer=0.2, Option_t * opt="");
  TStripDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X=0.,
           Int_t N_Strips=16, Double_t strip_width=0.2, Double_t inter_width=0.01, Double_t frame_width=0.2, Double_t dead_layer=0.2, Option_t * opt="");     
  ~TStripDetector();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area
  Int_t     GetStripFront(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
  Int_t     GetStripBack(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
  void      RotateX(Double_t);
  void      RotateY(Double_t);
  Double_t  GetThetaPixel(int stripf, int stripb);
  Double_t  GetPhiPixel(int stripf, int stripb);
  Double_t* GetThetaPhiPixel(int stripf, int stripb);
  Int_t     GetThetaPhiPixel(Double_t *,Double_t *, int stripf, int stripb);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TGraph*   GetGraphObject();
  TVector3  GetDetectorCenter(); // returns a TVector3 representing the center of the detector in the lab reference frame
  TVector3  GetPixelCenter(int stripf, int stripb); // returns a TVector3 representing the center of the pixel identified by a given strip front and a strip back in the lab reference frame
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
#ifdef GRAPHICAL_DEBUG
  void      ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
#endif
} ;

#endif
