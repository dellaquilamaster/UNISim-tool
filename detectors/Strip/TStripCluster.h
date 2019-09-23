/* *****************************************
 * 25/04/2019
 * This class handles a cluster of TStripDetector object
 * New draw functions added.
 * Created by: DELL'AQUILA DANIELE
 * Email:      daniele.dellaquila@irb.hr
 * *****************************************/
#ifndef TSTRIPCLUSTER_H
#define TSTRIPCLUSTER_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <TVector3.h>

#include "TStripDetector.h"

#define GRAPHICAL_DEBUG

class TStripCluster
{
public :
  TStripCluster();
  ~TStripCluster();
  
  void AddDetector(Double_t distance=15, Double_t theta_pos=0, Double_t phi_pos=0, Int_t N_Strips=16, 
		   Double_t strip_width=0.2, Double_t inter_width=0.01, Double_t frame_width=0.2, Double_t dead_layer=0.2, Option_t * opt="");  // Add a detector to the cluster facing the target
  void AddDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X=0.,
           Int_t N_Strips=16, Double_t strip_width=0.2, Double_t inter_width=0.01, Double_t frame_width=0.2, Double_t dead_layer=0.2, Option_t * opt="");  // Add a detector to the cluster with an arbitrary tilt
  void Clear();  // Clear the cluster
  bool CheckOverlap() const; // Check for detector overlap, if overlaps are found it returns true
  
  int Size() const;  // Get the number of telescopes in the cluster
  TStripDetector * GetDetector(int);  // Get a pointer to the TStripDetector at a certain index
  
  int IsInside(double theta, double phi, double x0=0, double y0=0, double z0=0) const;  // returns 1 if the particle is inside the cluster, 0 otherwise
  TStripDetector * GetDetector(double theta, double phi, double x0=0, double y0=0, double z0=0); // returns the pointer to the detector object that detects the particle, 0 if the particle is not inside the cluster
  int GetDetectorIndex(double theta, double phi, double x0=0, double y0=0, double z0=0) const; // returns the index of the detector object that detects the particle, -1 if the particle is not inside the cluster
  int GetDetectorStripFront(int numdet, double theta, double phi, double x0=0, double y0=0, double z0=0) const; // returns the strip front struck on a particolar detector, -1 if not in active area
  int GetDetectorStripBack(int numdet, double theta, double phi, double x0=0, double y0=0, double z0=0) const; // returns the strip back struck on a particolar detector, -1 if not in active area
  double GetDetectorPixelTheta(int numdet, int stripf, int stripb) const; // returns the theta angle (radiants) identified by a certain pixel of the detector numdet in the laboratory frame
  double GetDetectorPixelPhi(int numdet, int stripf, int stripb) const; // returns the phi angle (radiants) identified by a certain pixel of the detector numdet in the laboratory frame
  TVector3 GetDetectorPixelCenter(int numdet, int stripf, int stripb) const; // returns the vector identifying a certain pixel of the detector numdet with respect to the laboratory frame
  TVector3 GetDetectorCenter(int numdet) const; // returns the vector identifying the center of a given detector with respect to the origin    
  
  void Draw(Option_t * opt="", double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0) const; // Draw the whole cluster on the x-y plane
  void Draw3D(Option_t * opt="") const; // Draw the whole cluster on a 3D space
  
#ifdef GRAPHICAL_DEBUG
  void ShowImpactPoint(Double_t, Double_t, Double_t x0=0., Double_t y0=0., Double_t z0=0., double Xmin=0, double Xmax=0, double Ymin=0, double Ymax=0);
#endif
    
private :
  int fNumDetectors;
  std::vector <TStripDetector *> fTheDetectors;
    
};

#endif
