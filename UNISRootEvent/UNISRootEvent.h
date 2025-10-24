#ifndef UNISROOTEVENT
#define UNISROOTEVENT

#include <TROOT.h>

#include <UNISShared.h>

class UNISRootEvent
{
public :
  UNISRootEvent(); //!
  virtual ~UNISRootEvent(); //!
  
  Int_t fmulti;
  Int_t fmulti_detected;
  Int_t fsteps;
  Bool_t * fIsDetected; //[fmulti]
  Int_t * fnumdet; //[fmulti]
  Int_t * fnumpixel; //[fmulti]
  Double_t * fKinEnergy; //[fmulti]
  Double_t * fXDetHit; //[fmulti]
  Double_t * fYDetHit; //[fmulti]
  Double_t * fZDetHit; //[fmulti]
  Double_t * fKinEnergyAfterTarget; //[fmulti]
  Double_t * fKinEnergyOrigin; //[fmulti]
  Double_t * fThetaOrigin; //[fmulti]
  Double_t * fPhiOrigin; //[fmulti]
  Double_t * fThetaAfterTarget; //[fmulti]
  Double_t * fPhiAfterTarget; //[fmulti]
  Double_t * fZ; //[fmulti]
  Double_t * fA; //[fmulti]
  Double_t * fWeight; //[fsteps]

  ClassDef(UNISRootEvent,1);
};

#endif
