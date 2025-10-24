#include <UNISRootEvent.h>

//____________________________________________________
UNISRootEvent::UNISRootEvent() :
fmulti(0),
fsteps(0),
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
fThetaAfterTarget(new Double_t[gMaxEvtMulti]),
fPhiAfterTarget(new Double_t[gMaxEvtMulti]),
fZ(new Double_t[gMaxEvtMulti]),
fA(new Double_t[gMaxEvtMulti]),
fWeight(new Double_t[gMaxSteps])
{}

//____________________________________________________
UNISRootEvent::~UNISRootEvent()
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
  delete [] fThetaAfterTarget;
  delete [] fPhiAfterTarget;
  delete [] fZ;
  delete [] fA;
  delete [] fWeight;
}
