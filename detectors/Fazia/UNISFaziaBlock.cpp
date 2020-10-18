#include "UNISFaziaBlock.h"

//____________________________________________________
UNISFaziaBlock::UNISFaziaBlock(double theta_pos, double phi_pos, double displacement, Double_t pad_width, Double_t frame_width) :
fNumQuartets(4),
fNumQuartetRowColumn(sqrt(fNumQuartets)),
fQuartetWidth(pad_width*2+frame_width*4),
fQuartetHalfWidth(0.5*fQuartetWidth),
fQuartets(new UNISFaziaQuartet*[fNumQuartets]),
fNominalTheta(theta_pos),
fNominalPhi(phi_pos),
fNominalDistance(100),
fDisplacement(displacement)
{
  //
  //Creation of quartets
  for(int i=0; i<fNumQuartetRowColumn; i++) {
    for(int j=0; j<fNumQuartetRowColumn; j++) {      
      const double tilt_angle = atan(fQuartetHalfWidth/fNominalDistance);
      TVector3 AbsolutePosition(-fQuartetHalfWidth+2*i*fQuartetHalfWidth,fQuartetHalfWidth-2*j*fQuartetHalfWidth,fNominalDistance);
      fQuartets[i*fNumQuartetRowColumn+j] = new UNISFaziaQuartet(0.,0.,-100,pad_width,frame_width); //Quartet at 0 cm from the target
      fQuartets[i*fNumQuartetRowColumn+j]->RotateX(-tilt_angle+2*j*tilt_angle);
      fQuartets[i*fNumQuartetRowColumn+j]->RotateY(-tilt_angle+2*i*tilt_angle);
      fQuartets[i*fNumQuartetRowColumn+j]->Translate(AbsolutePosition.X(), AbsolutePosition.Y(), AbsolutePosition.Z());
      fQuartets[i*fNumQuartetRowColumn+j]->RotateZ(-phi_pos-180*TMath::DegToRad());
      fQuartets[i*fNumQuartetRowColumn+j]->RotateX(theta_pos);
      fQuartets[i*fNumQuartetRowColumn+j]->RotateZ(phi_pos+180*TMath::DegToRad());
      fQuartets[i*fNumQuartetRowColumn+j]->Translate(0., 0., fDisplacement);
    }
  }
  //
}

//____________________________________________________
UNISFaziaBlock::~UNISFaziaBlock()
{
  delete [] fQuartets; 
}

//
// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
//____________________________________________________
Int_t UNISFaziaBlock::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  for(Int_t i=0; i<fNumQuartets; i++)
  {
    if(fQuartets[i]->IsInside(theta_inc,phi_inc,x0,y0,z0)) return 1; 
  }
  
  //
  return 0;
  //
}

// 
// Returns an absolute number identifying the pixel fired according to the following scheme:
// back*num_pads + front
// If the particle is not inside the active area -> return value = -1
//____________________________________________________
Int_t UNISFaziaBlock::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  for(Int_t i=0; i<fNumQuartets; i++)
  {
    if(fQuartets[i]->IsInside(theta_inc,phi_inc,x0,y0,z0)) {
      const int pixel=fQuartets[i]->GetPixel(theta_inc,phi_inc,x0,y0,z0);
      //
      if(pixel<0) continue;
      //
      const int row = pixel/fNumQuartetRowColumn+(i/fNumQuartetRowColumn)*fNumQuartetRowColumn;
      const int column = pixel%fNumQuartetRowColumn+(i%fNumQuartetRowColumn)*fNumQuartetRowColumn;
      return row*fNumQuartetRowColumn*2+column;
      //
    } 
  }
  
  
  /*particle not inside the effective area*/
  return -1;
}

//
//draws the telescope on the X-Y plane
//____________________________________________________
void UNISFaziaBlock::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return; 
}

//____________________________________________________
TVector3 UNISFaziaBlock::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  for(Int_t i=0; i<fNumQuartets; i++)
  {
    if(fQuartets[i]->IsInside(theta_inc,phi_inc,x0,y0,z0)) return fQuartets[i]->GetImpactPointLab(theta_inc,phi_inc,x0,y0,z0); 
  }

  //
  return TVector3(0,0,0);
  //
}

//
//3D drawing function
//____________________________________________________
void UNISFaziaBlock::Draw3D(Option_t * draw_opt) const
{
  for(Int_t i=0; i<fNumQuartets; i++)
  {
    fQuartets[i]->Draw3D(draw_opt);
  }
}
