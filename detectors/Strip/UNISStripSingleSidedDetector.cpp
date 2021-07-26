#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "UNISStripSingleSidedDetector.h"

// standard constructor
UNISStripSingleSidedDetector::UNISStripSingleSidedDetector(Double_t distance, Double_t theta_pos, Double_t phi_pos, Int_t N_Strips, 
				   Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TTrueImpactPoint(0.,0.,0.),
TTelescopeImpactPoint(0.,0.,0.),
TTiltXAngle(-9999),
TTiltYAngle(-9999),
TStrips_number(N_Strips),                                                  
TStripTrue_width(strip_width),                                             
TStripTrue_semi(0.5*TStripTrue_width),                                     
TInter_width(inter_width),                                                 
TFrame_width(frame_width),                                                 
TDeadLayer(dead_layer),                                                    
TStripEffective_width(strip_width-inter_width),                            
TStripEffective_semi(0.5*TStripEffective_width),                           
TTelescopeEffective_semi(TStripTrue_semi*TStrips_number),                  
TTelescopeTrue_semi(TTelescopeEffective_semi+TDeadLayer+TFrame_width)
{
  //
  //Setting the position of the telescope center
  TTelescopeCenter.SetXYZ(0.,0.,distance);
  //Telescope's corners
  TTopLeftCorner.SetXYZ(TTelescopeTrue_semi,-TTelescopeTrue_semi,distance);
  TTopRightCorner.SetXYZ(TTelescopeTrue_semi,TTelescopeTrue_semi,distance);
  TBottomLeftCorner.SetXYZ(-TTelescopeTrue_semi,-TTelescopeTrue_semi,distance);
  TTopLeftXCorner= TTelescopeEffective_semi;
  TTopLeftYCorner= TTelescopeEffective_semi;
  //
  
  //
  //calculation of the X Y versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(0,0,1);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
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
UNISStripSingleSidedDetector::UNISStripSingleSidedDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y,
                   Int_t N_Strips, Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
UNISStripSingleSidedDetector(0., 0., 0., N_Strips, strip_width, inter_width, frame_width, dead_layer, opt)
{
  //
  //
  TTiltXAngle=tilt_X;
  TTiltYAngle=tilt_Y;
  //
  //
  //We operate now a series of rotations to reach the final position
  //First: Rotation about the X-axis of the titl_X angle
  RotateX(TTiltXAngle); 
  //Second: Rotation about the Y-axis of the titl_Y angle
  RotateY(TTiltYAngle);
  //Third: Translation to the desired position
  Translate(X0,Y0,Z0);
  //
}

void UNISStripSingleSidedDetector::RotateX(Double_t x_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateX(x_angle);
  //
  
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateX(x_angle); 
  TTopRightCorner  .RotateX(x_angle); 
  TBottomLeftCorner.RotateX(x_angle);
  //
  //calculation of the X Y versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateX(x_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  //
  Rotate3DX(x_angle);
  //  
  
  return;
}

void UNISStripSingleSidedDetector::RotateY(Double_t y_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateY(y_angle);
  //
  
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateY(y_angle); 
  TTopRightCorner  .RotateY(y_angle); 
  TBottomLeftCorner.RotateY(y_angle);
  //
  //calculation of the X Y versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateY(y_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  //
  Rotate3DY(y_angle);
  //  
  
  return;
}

void UNISStripSingleSidedDetector::RotateZ(Double_t z_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateZ(z_angle);
  //
  
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateZ(z_angle); 
  TTopRightCorner  .RotateZ(z_angle); 
  TBottomLeftCorner.RotateZ(z_angle);
  //
  //calculation of the X Y versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateZ(z_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  //
  Rotate3DZ(z_angle);
  //  
  
  return;
}


void UNISStripSingleSidedDetector::Translate(Double_t x, Double_t y, Double_t z)
{
  //
  TVector3 TranslationVector (x, y, z);
  
  //Translation of the telescope center
  TTelescopeCenter+=TranslationVector;
  
  //Translation of the corners
  TTopLeftCorner+=TranslationVector;
  TTopRightCorner+=TranslationVector;
  TBottomLeftCorner+=TranslationVector;
  
  //Calculation of detector reference frame vectors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag()); 
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  //
  Translate3D(x,y,z);
  //
  
  return;
}

// destructor
UNISStripSingleSidedDetector::~UNISStripSingleSidedDetector()
{}

// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
Int_t UNISStripSingleSidedDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  Double_t t; /*t parameter for the parametric equation of the particle direction*/
  Double_t l,m,n; /*particle direction parameters*/
  l=TMath::Sin(theta_inc)*TMath::Cos(phi_inc);
  m=TMath::Sin(theta_inc)*TMath::Sin(phi_inc);
  n=TMath::Cos(theta_inc);

  /*solution of the direction/plane interception*/
  if((t=-(Ta*x0+Tb*y0+Tc*z0+Td)/(Ta*l+Tb*m+Tc*n))<0) return 0;
  
  /*setting the true impact point*/
  TTrueImpactPoint.SetXYZ(x0+t*l,y0+t*m,z0+t*n);
    
  // setting he impact point in the telescope frame
  TTelescopeImpactPoint=TTrueImpactPoint-TTelescopeCenter;
  TImpactX=TTelescopeImpactPoint.Dot(TXversor);
  TImpactY=TTelescopeImpactPoint.Dot(TYversor); 
  
  /*inside condition*/
  if(TMath::Abs(TImpactX)<=TTelescopeEffective_semi && TMath::Abs(TImpactY)<=TTelescopeEffective_semi)
  {
    return 1;
  }
  return 0;
}

// 
// Returns an absolute number identifying the pixel fired according to the following scheme:
// back*num_strips + front
// If the particle is not inside the active area -> return value = -1
Int_t UNISStripSingleSidedDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  Int_t i,j;
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  /*If not the particle is inside the telescope*/

  // Impact point coordinates with respect to the Top Left corner
  // These quantities are always positive within the surface of the detector
  TImpactXprime=TTopLeftXCorner-TImpactX;
  TImpactYprime=TImpactY+TTopLeftYCorner;
  j=Int_t(TImpactXprime/TStripTrue_width);
    
  /*check if the particle is inside the effective area*/
  if(std::fabs(TImpactXprime-(2*j+1)*TStripTrue_semi)<=TStripEffective_semi) 
  {
    return j;
  }
  /*particle not inside the effective area*/
  return -1;
}

//draws the telescope on the X-Y plane
void UNISStripSingleSidedDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return; 
}

//3D drawing function
void UNISStripSingleSidedDetector::Draw3D(Option_t * draw_opt) const
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

void UNISStripSingleSidedDetector::Generate3D()
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
  fDetector = new TGeoVolumeAssembly("DetectorQuartet");
  fDetectorMatrix = new TGeoHMatrix("DetectorTransformationMatrix");
  //
   
  //
  fFrame       = new TGeoVolume("frame_volume",new TGeoBBox(TTelescopeTrue_semi, TTelescopeTrue_semi+0.5, 0.05));
  fStrip       = new TGeoVolume("strip_volume",new TGeoBBox(TStripEffective_semi, TStripTrue_semi*TStrips_number, 0.1));
  fStripGround = new TGeoVolume("strip_ground_volume",new TGeoBBox(TStripTrue_semi*TStrips_number, TStripTrue_semi*TStrips_number, 0.05));
  fFrame->SetLineColor(kAzure+1);
  fStrip->SetLineColor(kGray);
  fStripGround->SetLineColor(kGray);
  //
  
  //
  //Adding frame to mother volume
  fDetector->AddNode(fFrame,0,new TGeoTranslation(0., -0.5, 0.));
  fDetector->AddNode(fStripGround,0,new TGeoTranslation(0., 0., 0.05));
  //
  
  //
  //Pixels
  for(Int_t i=0; i<TStrips_number; i++)
  {  
    fDetector->AddNode(fStrip,i,new TGeoTranslation(-(TStrips_number-(2*i+1))*TStripTrue_semi,0,0));
  }
  //
}  

void UNISStripSingleSidedDetector::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(x_angle*TMath::RadToDeg());
}

void UNISStripSingleSidedDetector::Rotate3DY(Double_t y_angle)
{
  fDetectorMatrix->RotateY(y_angle*TMath::RadToDeg());
}

void UNISStripSingleSidedDetector::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISStripSingleSidedDetector::Translate3D(Double_t x, Double_t y, Double_t z)
{  
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(x,y,z));
}

// returns the TVector3 of the telescope center
TVector3 UNISStripSingleSidedDetector::GetDetectorCenter()
{
  return TTelescopeCenter; 
}

TVector3 UNISStripSingleSidedDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TImpactX*TXversor+TImpactY*TYversor+TTelescopeCenter;
}

#ifdef GRAPHICAL_DEBUG
void UNISStripSingleSidedDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{}
#endif
