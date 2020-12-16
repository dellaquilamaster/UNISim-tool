#include <UNISDetectedRootEvent.h>

//____________________________________________________
UNISDetectedRootEvent::UNISDetectedRootEvent() :
fmulti(0),
fIsDetected(new Bool_t[gMaxEvtMulti]),
fnumdet(new Int_t[gMaxEvtMulti]),
fnumpixel(new Int_t[gMaxEvtMulti]),
fKinEnergy(new Double_t[gMaxEvtMulti]),
fXDetHit(new Double_t[gMaxEvtMulti]),
fYDetHit(new Double_t[gMaxEvtMulti]),
fZDetHit(new Double_t[gMaxEvtMulti]),
fKinEnergyAfterTarget(new Double_t[gMaxEvtMulti]),
fKinEnergyOrigin(new Double_t[gMaxEvtMulti]),
fThetaOrigin(new Double_t[gMaxEvtMulti]),
fPhiOrigin(new Double_t[gMaxEvtMulti]),
fThetaDetected(new Double_t[gMaxEvtMulti]),
fPhiDetected(new Double_t[gMaxEvtMulti]),
fKinEnergyOriginCms(new Double_t[gMaxEvtMulti]),
fThetaOriginCms(new Double_t[gMaxEvtMulti]),
fZ(new Int_t[gMaxEvtMulti]),
fA(new Int_t[gMaxEvtMulti])
{}

//____________________________________________________
UNISDetectedRootEvent::~UNISDetectedRootEvent()
{
  delete [] fIsDetected;
  delete [] fnumdet;
  delete [] fnumpixel;
  delete [] fKinEnergy;
  delete [] fXDetHit;
  delete [] fYDetHit;
  delete [] fZDetHit;
  delete [] fKinEnergyAfterTarget;
  delete [] fKinEnergyOrigin;
  delete [] fThetaOrigin;
  delete [] fPhiOrigin;
  delete [] fThetaDetected;
  delete [] fPhiDetected;
  delete [] fKinEnergyOriginCms;
  delete [] fThetaOriginCms;
  delete [] fZ;
  delete [] fA;
}
