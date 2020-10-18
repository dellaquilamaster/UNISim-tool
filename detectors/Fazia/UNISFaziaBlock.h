/* *****************************************
 * 25/05/2020
 * Class to handle a simple Fazia block
 * facing the target perpendicularly at a distance 
 * of 100 cm.
 * Each quartet in the block is facing the target perpendicularly at 100 cm distance.
 * A block is a grouping of 4 quartets.
 * The class is derived from TStrip
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef UNISFAZIABLOCK_H
#define UNISFAZIABLOCK_H

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

#include "../DetectionSetup/UNISDetectionUnit.h"
#include "UNISFaziaQuartet.h"

#define GRAPHICAL_DEBUG

/* *****************************************
 * 
 * Example
 * 
 * 4 quartets each containing 2 x 2 pads

 *  __________ __________
 * | 0  | 1  || 0  | 1  |                               
 * |----------|----------                               
 * | 2  | 3  || 2  | 3  |                                       
 * |---------||---------|                                       
 * | 0  | 1  || 0  | 1  |                                                  
 * |----------|----------
 * | 2  | 3  || 2  | 3  |
 * |---------||---------|
 *
 * Each quartet is facing the beam perpendicularly at a distance of 100 cm
 *
 * NOTE: The absolute number of pixel is considered from top left to bottom right.
 *
 * 
 * The reference frame used for the calculations has the z as the beam axis. In the upstream view, 
 * x is vertical and y goes towards the right hand side.
 * By defayult, detectors at phi=0 are located on the horizontal plane at the right side of the beam (upstream)
 * *****************************************/

class UNISFaziaBlock : public UNISDetectionUnit
{
public :
  UNISFaziaBlock(double theta_pos, double phi_pos, double displacement=0., Double_t pad_width=2., Double_t frame_width=0.1);
  ~UNISFaziaBlock();
  
  Int_t     IsInside(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns 1 if the particle is inside the detector
  Int_t     GetPixel(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; //returns an absolute number identifying the number of pad fired, -1 if not inside the active area

  void      Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const override;
  void      Draw3D(Option_t * opt="") const override;
  
  TVector3  GetImpactPointLab(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0.) override; // Get the impact point in the lab reference frame
    
private :
  const int fNumQuartets;
  const int fNumQuartetRowColumn;
  const double fQuartetWidth;
  const double fQuartetHalfWidth;
  UNISFaziaQuartet ** fQuartets;
  double fNominalTheta;
  double fNominalPhi;
  double fNominalDistance;
  double fDisplacement;
  
};

#endif
