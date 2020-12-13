#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "UNISCollimatedSilicon.h"

/*
 * In the lab frame, looking from the target (downstream)
 * 
 *    ^ y-axis
 *    |
 *    |
 * <---
 * x-axis
 * 
 * y-axis penetrates the screen.
 * 
 * In the pad frame
 * 
 *   ^ y'-axis
 *   |
 *   |
 *   ---> x'-axis
 * 
 */

// standard constructor
UNISCollimatedSilicon::UNISCollimatedSilicon(Double_t distance, Double_t theta_pos, Double_t phi_pos, Double_t collimator_inner_radius, Double_t collimator_outer_radius) :
fXlabversor(1,0,0),
fYlabversor(0,1,0),
fZlabversor(0,0,1),
fCenter(0,0,0),
fLabImpactPoint(0.,0.,0.),
fPadImpactPoint(0.,0.,0.),
fTiltXAngle(-9999),
fTiltYAngle(-9999),
fTiltZAngle(-9999),                                           
fCollimatorInnerRadius(collimator_inner_radius),                                                                                       
fCollimatorOuterRadius(collimator_outer_radius),                            
fImpactX(0),
fImpactY(0)
{
  //
  //Setting the position of the telescope center
  fCenter.SetXYZ(0.,0.,distance);
  //Telescope's corners
  fTopLeftReference.SetXYZ(fCollimatorInnerRadius,fCollimatorInnerRadius,distance);
  fTopRightReference.SetXYZ(-fCollimatorInnerRadius,fCollimatorInnerRadius,distance);
  fBottomLeftReference.SetXYZ(fCollimatorInnerRadius,-fCollimatorInnerRadius,distance);
  fBottomRightReference.SetXYZ(-fCollimatorInnerRadius,-fCollimatorInnerRadius,distance);
  //
  
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftReference-fBottomLeftReference);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightReference-fTopLeftReference);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(0,0,1);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  Generate3D();
  //
  
  //
  Translate3D(0,0,distance);
  //
    
  //
  //We operate now a series of rotations to reach the final position
  //First: Rotation about the Z-axis of a quantity (-phi)
  RotateZ(-phi_pos-180*TMath::DegToRad()); 
  //Second: Rotation about the X-axis of a quantity (theta)
  RotateX(theta_pos);
  //Third: Rotation about the Z-axis of a quantity (phi)
  RotateZ(phi_pos+180*TMath::DegToRad());
  //  
}  

// constructor with arbitrary tilt angle, center position (X0, Y0, Z0)
// tilt_X = tilt angle with respect to the X-axis (vertical)
UNISCollimatedSilicon::UNISCollimatedSilicon(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y, Double_t tilt_Z, Double_t collimator_inner_radius, Double_t collimator_outer_radius) :
UNISCollimatedSilicon(0., 0., 0.)
{
  //
  //
  fTiltXAngle=tilt_X;
  fTiltYAngle=tilt_Y;
  fTiltZAngle=tilt_Z;
  //
  //
  //We operate now a series of rotations to reach the final position
  //First: Rotation about the X-axis of the titl_X angle
  RotateX(fTiltXAngle); 
  //Second: Rotation about the Y-axis of the titl_Y angle
  RotateY(fTiltYAngle);
  //Third: Rotation about the Z-axis of the titl_Z angle
  RotateZ(fTiltZAngle);
  //Fourth: Translation to the desired position
  Translate(X0,Y0,Z0);
  //
}

void UNISCollimatedSilicon::RotateX(Double_t x_angle)
{
  /*Rotation of the telescope center to the input position*/
  fCenter.RotateX(x_angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftReference   .RotateX(x_angle); 
  fTopRightReference  .RotateX(x_angle); 
  fBottomLeftReference.RotateX(x_angle);
  fBottomRightReference.RotateX(x_angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftReference-fBottomLeftReference);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightReference-fTopLeftReference);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateX(x_angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  Rotate3DX(x_angle);
  //  
  
  return;
}

void UNISCollimatedSilicon::RotateY(Double_t y_angle)
{
  /*Rotation of the telescope center to the input position*/
  fCenter.RotateY(y_angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftReference   .RotateY(y_angle); 
  fTopRightReference  .RotateY(y_angle); 
  fBottomLeftReference.RotateY(y_angle);
  fBottomRightReference.RotateY(y_angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftReference-fBottomLeftReference);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightReference-fTopLeftReference);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateY(y_angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  Rotate3DY(y_angle);
  //  
  
  return;
}

void UNISCollimatedSilicon::RotateZ(Double_t z_angle)
{
  /*Rotation of the telescope center to the input position*/
  fCenter.RotateZ(z_angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftReference   .RotateZ(z_angle); 
  fTopRightReference  .RotateZ(z_angle); 
  fBottomLeftReference.RotateZ(z_angle);
  fBottomRightReference.RotateZ(z_angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftReference-fBottomLeftReference);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightReference-fTopLeftReference);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateZ(z_angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  Rotate3DZ(z_angle);
  //  
  
  return;
}


void UNISCollimatedSilicon::Translate(Double_t x, Double_t y, Double_t z)
{
  //
  TVector3 TranslationVector (x, y, z);
  
  //Translation of the telescope center
  fCenter+=TranslationVector;
  
  //Translation of the corners
  fTopLeftReference+=TranslationVector;
  fTopRightReference+=TranslationVector;
  fBottomLeftReference+=TranslationVector;
  fBottomRightReference+=TranslationVector;
  
  //calculation of the X Y versors
  fYversor=(fTopLeftReference-fBottomLeftReference);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightReference-fTopLeftReference);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  Translate3D(x,y,z);
  //
  
  return;
}

// destructor
UNISCollimatedSilicon::~UNISCollimatedSilicon()
{}

//
// returns 1 if the particle is inside the photodiode, 0 if not.
// this function sets also the impact point XY coordinates
Int_t UNISCollimatedSilicon::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  Double_t t; /*t parameter for the parametric equation of the particle direction*/
  Double_t l,m,n; /*particle direction parameters*/
  l=TMath::Sin(theta_inc)*TMath::Cos(phi_inc);
  m=TMath::Sin(theta_inc)*TMath::Sin(phi_inc);
  n=TMath::Cos(theta_inc);

  /*solution of the direction/plane interception*/
  if((t=-(fa*x0+fb*y0+fc*z0+fd)/(fa*l+fb*m+fc*n))<0) return 0;
  
  /*setting the true impact point*/
  fLabImpactPoint.SetXYZ(x0+t*l,y0+t*m,z0+t*n);
    
  // setting he impact point in the telescope frame
  fPadImpactPoint=fLabImpactPoint-fCenter;
  fImpactX=fPadImpactPoint.Dot(fXversor);
  fImpactY=fPadImpactPoint.Dot(fYversor); 
  
  /*inside condition*/
  if(sqrt(fImpactX*fImpactX+fImpactY*fImpactY)<=fCollimatorOuterRadius)
  {
    return 1;
  }
  return 0;
}

// 
// Returns 0 if the particle is inside the active area:
// If the particle is not inside the active area -> return value = -1
Int_t UNISCollimatedSilicon::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  /*If not the particle is inside the telescope*/
  
  /*check if the particle is inside the effective area*/
  if(sqrt(fImpactX*fImpactX+fImpactY*fImpactY)<=fCollimatorInnerRadius) 
  {
    return 0;
  }
  //particle not inside the effective area
  return -1;
}

//draws the telescope on the X-Y plane
void UNISCollimatedSilicon::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return; 
}

//3D drawing function
void UNISCollimatedSilicon::Draw3D(Option_t * draw_opt) const
{
  //
  if(strstr(draw_opt,"SAME")==0 && strstr(draw_opt,"same")==0) {
    // 
    TGeoTube * Z_axis_line_shape = new TGeoTube("Z_axis_line_shape",0, 0.2, 15);
    TGeoCone * Z_axis_end_shape = new TGeoCone("Z_axis_end_shape",2,0,1,0,0);
    TGeoTranslation * Z_axis_arrow_trans = new TGeoTranslation("Z_axis_arrow_trans",0,0,15);
    Z_axis_arrow_trans->RegisterYourself();
    TGeoCompositeShape * Z_axis_shape = new TGeoCompositeShape("","(Z_axis_line_shape+Z_axis_end_shape:Z_axis_arrow_trans)");
    TGeoVolume * Z_axis_volume = new TGeoVolume("Z_axis_volume", Z_axis_shape);
    TGeoHMatrix * Z_axis_matrix = new TGeoHMatrix();
    Z_axis_volume->SetLineColor(kBlue);
    gGeoManager->GetMasterVolume()->AddNode(Z_axis_volume,0,Z_axis_matrix);
    Z_axis_matrix->Multiply(TGeoTranslation(0,0,15));
    
    TGeoTube * X_axis_line_shape = new TGeoTube("X_axis_line_shape",0, 0.2, 15);
    TGeoCone * X_axis_end_shape = new TGeoCone("X_axis_end_shape",2,0,1,0,0);
    TGeoTranslation * X_axis_arrow_trans = new TGeoTranslation("X_axis_arrow_trans",0,0,15);
    X_axis_arrow_trans->RegisterYourself();
    TGeoCompositeShape * X_axis_shape = new TGeoCompositeShape("(X_axis_line_shape+X_axis_end_shape:X_axis_arrow_trans)");
    TGeoVolume * X_axis_volume = new TGeoVolume("X_axis_volume", X_axis_shape);
    TGeoHMatrix * X_axis_matrix = new TGeoHMatrix();
    X_axis_volume->SetLineColor(kRed);
    gGeoManager->GetMasterVolume()->AddNode(X_axis_volume,0,X_axis_matrix);
    X_axis_matrix->RotateY(90);
    X_axis_matrix->Multiply(TGeoTranslation(0,0,15));
    
    TGeoTube * Y_axis_line_shape = new TGeoTube("Y_axis_line_shape",0, 0.2, 15);
    TGeoCone * Y_axis_end_shape = new TGeoCone("Y_axis_end_shape",2,0,1,0,0);
    TGeoTranslation * Y_axis_arrow_trans = new TGeoTranslation("Y_axis_arrow_trans",0,0,15);
    Y_axis_arrow_trans->RegisterYourself();
    TGeoCompositeShape * Y_axis_shape = new TGeoCompositeShape("(Y_axis_line_shape+Y_axis_end_shape:Y_axis_arrow_trans)");
    TGeoVolume * Y_axis_volume = new TGeoVolume("Y_axis_volume", Y_axis_shape);
    TGeoHMatrix * Y_axis_matrix = new TGeoHMatrix();
    Y_axis_volume->SetLineColor(kGreen+1);
    gGeoManager->GetMasterVolume()->AddNode(Y_axis_volume,0,Y_axis_matrix);
    Y_axis_matrix->RotateX(-90);
    Y_axis_matrix->Multiply(TGeoTranslation(0,0,15));
    //
  }
  //
                     
  //
  gGeoManager->GetMasterVolume()->AddNode(fDetector,1,fDetectorMatrix);
  //
  
  //
  std::string option_string(draw_opt);
  gGeoManager->GetMasterVolume()->Draw(option_string.find("ogl")!=std::string::npos ? "ogl" : "");
  //
}

void UNISCollimatedSilicon::Generate3D()
{
  //
  if(!gGeoManager) {
    new TGeoManager();
    TGeoVolume *TheMotherVolume = gGeoManager->MakeBox("TheMotherVolume",0,50,50,50);
    TheMotherVolume->SetLineColor(kBlack);
    gGeoManager->SetTopVisible(kFALSE); // the TOP is invisible
    gGeoManager->SetTopVolume(TheMotherVolume);
  }
  
  //
  fDetector = new TGeoVolumeAssembly("Detector");
  fDetectorMatrix = new TGeoHMatrix("DetectorTransformationMatrix");
  //
   
  //
  fCollimator   = new TGeoVolume("collimator_volume",new TGeoTube(fCollimatorInnerRadius, fCollimatorOuterRadius, 0.2));
  fActiveArea   = new TGeoVolume("active_area_volume",new TGeoTube(0., fCollimatorInnerRadius, 0.1));
  fCollimator->SetLineColor(kYellow+2);
  fActiveArea->SetLineColor(kGray);
  //
  
  //
  //Adding to mother volume
  fDetector->AddNode(fCollimator,0,new TGeoTranslation(0., 0., 0.));
  fDetector->AddNode(fActiveArea,0,new TGeoTranslation(0., 0., 0.));
  //
}  

void UNISCollimatedSilicon::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(x_angle*TMath::RadToDeg());
}

void UNISCollimatedSilicon::Rotate3DY(Double_t y_angle)
{
  fDetectorMatrix->RotateY(y_angle*TMath::RadToDeg());
}

void UNISCollimatedSilicon::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISCollimatedSilicon::Translate3D(Double_t x, Double_t y, Double_t z)
{  
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(x,y,z));
}

// returns the TVector3 of the telescope center
TVector3 UNISCollimatedSilicon::GetDetectorCenter()
{
  return fCenter; 
}

TVector3 UNISCollimatedSilicon::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return fLabImpactPoint;
}
