/* *****************************************
 * TLampDetector
 * 22/12/2015
 * C++ Object made to simulate a 
 * wedge of a lamp-like detector.
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef _TLAMPWEDGEDETECTOR
#define _TLAMPWEDGEDETECTOR

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
#include <TEveGeoShape.h>
#include <TEveTrans.h>
#include <TEveText.h>
#include <TGLViewer.h>

#include "../DetectionSetup/TDetectionUnit.h"

#define GRAPHICAL_DEBUG

/* *****************************************
 * The strips are annular region going from 0 (close to the beam line) to N_Strips-1 (outer strip).
 * The detector region has the shape described below.
 * 
 *           
 *          ^ x'
 *__________|___________    
 *\    _____|______    / 
 * \  /     |      \  /     
 *  \/      |silicon\/
 *   \      |area   /            
 *    \     |      /             
 *     \    |     /              
 *      \___|____/   
 *       \__|__ /
 *          |
 *          |
 *          |
 *          o The reference (beam axis) -------> y'
 * 
 * For the calculation we use the polar coordinates.
 * x' = rho*cos(theta)
 * y' = rho*sin(theta)
 * where theta is the angle formed between the vector 
 * identifying an impact point on the detector with respect 
 * to the ideal center and the vertical x' axis (lying on the detector and 
 * also symmetry axis of the detector).
 * y' is the horizontal axis in the reference frame of the detector.
 * (The latter is tangent to each strip at y'=0 and lies on the plane of the detector.
 * 
 * *****************************************/

class TLampWedgeDetector : public TDetectionUnit
{
private:
  TVector3   TXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   TYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   TZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3   TDetectorReference; //The bottom center point of the detector
  TVector3   TDetectorNominalReference; //The nominal reference of the detector, used for the calculation. This point always lies on the beam line.
  TVector3   TDetectorTopReference;  // A reference point to identify the vertical x' axis
  TVector3   TDetectorRightReference;  // A reference point to identify the horizontal y' axis
  TVector3   TXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   TYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   TTrueImpactPoint; /*True impact point coordinates in the laboratory frame*/
  TVector3   TDetectorImpactPoint; /*Impact point coordinates in the detector frame*/
  Double_t   TDistanceBeamAxis; //Detector position along the beam axis
  Double_t   TAzimuthalAngle; //Phi angle where the frame is placed
  Double_t   TTiltAngle; /*tilt angle with respect to the vertical (X) axis.*/
  Int_t      TStrips_number; /*Total number of strips*/
  Double_t   TStripTrue_width; /*width of each strip in mm*/ 
  Double_t   TStripTrue_semi; /*half-width of the strip including the interstrip*/
  Double_t   TInter_width; /*Width of the inter-strip of each strip*/
  Double_t   TStripEffective_width; /*width of the strip excluding the interstrip*/
  Double_t   TStripEffective_semi; /*half-width of the strip excluding the interstrip*/
  Double_t   TInnerNominal_radius; //Radius of the inner side of the most inner strip when the detector is placed in its nominal position
  Double_t   TStripNominalCoverageAngle; //Nominal coverage angle of the strips.
  Double_t * TStripCoverageAngle; //Strip-by-strip array of maximum angle (with respect to the x' axis in the primed reference frame) seen by a strip
  Double_t   TFrameCoverageAngle; //Maximum angle covered by the frame
  Double_t   TBottomFrame; //Size of the bottom frame, this is also placed on the top
  Double_t   TBottomFrame_distance; //Distance from the bottom frame to the beam line
  Double_t   TTopFrame_distance; //Distance from the top frame to the beam line
  Double_t   TNominalDistanceBeamLine; //Nominal distance from the bottom edge of the frame to the beam line
  Double_t   TNominalDistanceTopBeamLine; //Nominal distance from the top edge of the frame to the beam line
  Double_t * TStripRadius; //Radius of the i-th strip
  Double_t   TImpactX; /*X coordinate of the impact point on the telescope surface with respect to the detector reference in the telescope frame*/
  Double_t   TImpactY; /*Y coordinate of the impact point on the telescope surface with respect to the detector reference in the telescope frame*/
  Double_t   Ta; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   Tb;
  Double_t   Tc;
  Double_t   Td;
  
public:
  TLampWedgeDetector(Double_t distance=15, Double_t phi_pos=0, Double_t tilt=0., Double_t bottom_frame_distance=2, Int_t N_Strips=16, Double_t strip_width=0.5, Double_t inter_width=0.01, Double_t nominal_radius=5., 
                     Double_t nominal_coverage=42., Double_t nominal_frame_coverage=45, Double_t bottom_frame=1., Option_t * opt="");    
  ~TLampWedgeDetector();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area;

  void      RotateX(Double_t); //Rotation of the whole detector around the Y-axis 
  void      RotateY(Double_t); //Rotation of the whole detector around the X-axis
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TGraph*   GetGraphObject();
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
#ifdef GRAPHICAL_DEBUG
  void      ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
#endif
} ;

#endif
