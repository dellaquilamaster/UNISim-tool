/* *****************************************
 * UNISS3Detector
 * 06/04/2023
 * C++ Object made to simulate a 
 * Micron Semiconductors S3 detector
 * with central hole.
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@unina.it
 * *****************************************/
#ifndef UNISS3DETECTOR_H
#define UNISS3DETECTOR_H

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
#include <TGeoTube.h>
#include <TGeoTrd2.h>
#include <TGeoCompositeShape.h>

#include "../DetectionSetup/UNISDetectionUnit.h"

#define GRAPHICAL_DEBUG

/* *****************************************
 * The strips are annular regions going from 0 (close to the beam line) to TAnnularStrips_number-1 (outer strip).
 * The detector has a central hole and is meant to be placed looking at the target 
 * at a given theta and phi position and a given distance from the target.
 * From the center of the central hole, the beginning of the silicon region 
 * starts at TInnerRadius and ends at TOuterRadius. The active area starts at 
 * TInnerActiveAreaRadius and ends at TOuterActiveAreaRadius. There are TAnnularStrips_number=24 annular strips 
 * and TRadialStrips_number=32 radial strips.
 * 
 * Annular strips have a pitch of TAnnularStripTrue_width=1 mm, and an effective size of TAnnularStripEffective_width=0.886 mm. 
 * Radial strips have a pitch of TRadialStripCoverageAngle=11.25 deg.
 * 
 * The following picture represents a sector of the detector. The actual detector 
 * is a CD-rom around the central axis.
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
 * Front side is segmented into annular strips of minimum radius (inner border) TInnerRadius and maximum radius (outer border) 
 * TOuterRadius. Back side is sgemented into radial strips.
 * 
 * *****************************************/

class UNISS3Detector : public UNISDetectionUnit
{
private:
  TVector3   TXlabversor; // X-axis versor in the lab frame (vertical axis)
  TVector3   TYlabversor; // Y-axis versor in the lab frame (horizontal axis)
  TVector3   TZlabversor; // Z-axis versor in the lab frame (beam axis)
  TVector3   TDetectorReference; //The center of the detector
  TVector3   TDetectorTopReference;  // A reference point to identify the vertical x' axis
  TVector3   TDetectorRightReference;  // A reference point to identify the horizontal y' axis
  TVector3   TXversor; /*X-axis versor in the telescope frame (X' axis)*/
  TVector3   TYversor; /*Y-axis versor in the telescope frame (Y' axis)*/
  TVector3   TTrueImpactPoint; /*True impact point coordinates in the laboratory frame*/
  TVector3   TDetectorImpactPoint; /*Impact point coordinates in the detector frame*/
  Double_t   TDistance; //Distance from detector center to target
  Double_t   TPolarAngle; //Theta angle where the frame is placed (geometrical center of the detector)
  Double_t   TAzimuthalAngle; //Phi angle where the frame is placed (geometrical center of the detector)
  Int_t      TAnnularStrips_number; //Total number of annular strips (24)
  Int_t      TRadialStrips_number; //Total number of radial strips (32)
  Double_t   TAnnularStripTrue_width; //Width of each annular strip (1.0 mm)
  Double_t   TAnnularStripTrue_semi; //Half-width of each annular strip
  Double_t   TInterWidth; //Width of the inter-strip
  Double_t   TAnnularStripEffective_width; //Width of the effective area of each annular strip (0.886 mm)
  Double_t   TAnnularStripEffective_semi;  //Half-width of the effective area of each annular strip
  Double_t   TInnerRadius; //Inner radius of the silicon region (20 mm)
  Double_t   TOuterRadius; //Outer radius of the silicon region (76 mm)
  Double_t   TInnerActiveAreaRadius; //Inner radius of the effective area of the silicon region (22 mm)
  Double_t   TOuterActiveAreaRadius; //Outer radius of the effective area of the silicon region (70 mm)
  Double_t   TRadialStripCoverageAngle; //Nominal coverage angle of the radial strips.
  Double_t   TFrameWidth; //Width of the frame
  Double_t * TRadialStripMinimumEffectiveAngle; //Strip-by-strip array of minimum angle (with respect to the x' axis in the primed reference frame) seen by a radial strip excluding interstrip region
  Double_t * TRadialStripMaximumEffectiveAngle; //Strip-by-strip array of maximum angle (with respect to the x' axis in the primed reference frame) seen by a radial strip excluding interstrip region
  Double_t * TAnnularStripRadius; //Radius of the i-th annular strip  
  Double_t   TImpactX; /*X coordinate of the impact point on the telescope surface with respect to the detector reference in the telescope frame*/
  Double_t   TImpactY; /*Y coordinate of the impact point on the telescope surface with respect to the detector reference in the telescope frame*/
  Double_t   Ta; /*parameters for the equation idetifying the plane of the detector*/
  Double_t   Tb;
  Double_t   Tc;
  Double_t   Td;
  
public:
  UNISS3Detector(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0., Double_t inter_width=0.0114);
  ~UNISS3Detector();

  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pixel fired, -1 if not inside the active area;

  void      RotateX(Double_t); //Rotation of the whole detector around the X-axis 
  void      RotateZ(Double_t); //Rotation of the whole detector around the Z-axis
  void      TranslateLongitudinal(Double_t);
  void      Translate(Double_t, TVector3);
  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  TGraph*   GetGraphObject();
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
  
private :
  TGeoVolume * fFrame;
  TGeoVolume * fSilicon;
  TGeoVolume ** fStripAnnular;
  TGeoVolume ** fStripRadial;

  //
  void Generate3D(Double_t, Double_t);
  void Rotate3DX(Double_t);
  void Rotate3DZ(Double_t);
  void TranslateLongitudinal3D(Double_t);
  void Translate3D(Double_t, TVector3);
  
#ifdef GRAPHICAL_DEBUG
  void      ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.);
#endif
} ;

#endif
