#include "UNISLampWedgeDetector.h"


//____________________________________________________
UNISLampWedgeDetector::UNISLampWedgeDetector(Double_t distance, Double_t phi_pos, Double_t tilt, Double_t bottom_frame_distance, Int_t N_Strips, Double_t strip_width, Double_t inter_width, Double_t nominal_radius, 
                                             Double_t nominal_coverage, Double_t nominal_frame_coverage, Double_t bottom_frame, Option_t * opt) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TDistanceBeamAxis(distance),
TAzimuthalAngle(phi_pos),
TTiltAngle(tilt),
TBottomFrame_distance(bottom_frame_distance),
TStrips_number(N_Strips),                                                  
TStripTrue_width(strip_width),                                             
TStripTrue_semi(0.5*TStripTrue_width),                                     
TInter_width(inter_width),   
TStripEffective_width(TStripTrue_width-TInter_width),
TStripEffective_semi(TStripEffective_width/2.),
TInnerNominal_radius(nominal_radius),
TStripNominalCoverageAngle(nominal_coverage/2.*TMath::DegToRad()),
TStripCoverageAngle(new Double_t[TStrips_number]),
TFrameCoverageAngle(nominal_frame_coverage/2.*TMath::DegToRad()),                                                    
TBottomFrame(bottom_frame),
TTopFrame_distance(TBottomFrame_distance+TBottomFrame+TStripTrue_width*TStrips_number+TBottomFrame),
TNominalDistanceBeamLine(TInnerNominal_radius-TBottomFrame),
TNominalDistanceTopBeamLine(TNominalDistanceBeamLine+2*TBottomFrame+TStrips_number*TStripTrue_width),
TStripRadius(new Double_t[TStrips_number])
{
  //
  //Placing the detector plane orthogonal to the beam axis and with is bottom center touching the beam line
  //The plane is identified by the three reference points (origin on the beam-line and 2 additional points)
  TDetectorReference.SetXYZ(0.,0.,0.);
  TDetectorTopReference=TDetectorReference+TYlabversor;
  TDetectorRightReference=TDetectorReference-TXlabversor;
  //
  //Moving reference point according to the bottom frame displacement
  TDetectorReference+=TBottomFrame_distance*TYlabversor;
  TDetectorTopReference+=TBottomFrame_distance*TYlabversor; 
  TDetectorRightReference+=TBottomFrame_distance*TYlabversor;
  //
  
  //
  //Generating strips
  for(int i=0; i<TStrips_number; i++) {
    TStripRadius[i]=TInnerNominal_radius+(2*i+1)*TStripTrue_semi;
    TStripCoverageAngle[i]=TStripNominalCoverageAngle;
    //
    //Last 3 strips have a narrower coverage angle respectively 36.4 deg (86.67%), 29.2 deg (69.52%) and 18.9 deg (45%)
    if(i==TStrips_number-3) TStripCoverageAngle[i]*=0.8667;
    if(i==TStrips_number-2) TStripCoverageAngle[i]*=0.6952;
    if(i==TStrips_number-1) TStripCoverageAngle[i]*=0.4500;
    //
  }
  //
  
  //
  //calculation of the X Y versors after rotations
  TXversor=(TDetectorTopReference-TDetectorReference);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TDetectorRightReference-TDetectorReference);
  TYversor*=(1./TYversor.Mag());
  //
  
  //
  //Calculation of the reference origin
  TDetectorNominalReference=TDetectorReference-TNominalDistanceBeamLine*TXversor;
  //
    
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(0,0,1);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TDetectorReference.X()+Tb*TDetectorReference.Y()+Tc*TDetectorReference.Z());
  //
  
  //
  Generate3D(0.,0.);
  //
  
  //
  //Rotation X-tilt
  RotateX(tilt);
  //
  
  //
  //Translation to the desired distance
  TranslateLongitudinal(TDistanceBeamAxis);
  //
  
  //
  //Azimuthal rotation
  RotateZ(phi_pos+180*TMath::DegToRad());
  //
}

//____________________________________________________
UNISLampWedgeDetector::~UNISLampWedgeDetector()
{}

//____________________________________________________
void UNISLampWedgeDetector::RotateX(Double_t angle)
{
  //
  //Rotation of the telescope reference frame to the input position
  TDetectorReference.RotateX(angle);
  TDetectorTopReference.RotateX(angle);
  TDetectorRightReference.RotateX(angle);
  //
  
  //
  //calculation of the X Y versors
  TXversor=(TDetectorTopReference-TDetectorReference);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TDetectorRightReference-TDetectorReference);
  TYversor*=(1./TYversor.Mag());    
  //
  
  //
  //Calculation of the reference origin
  TDetectorNominalReference=TDetectorReference-TNominalDistanceBeamLine*TXversor;
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateX(angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TDetectorReference.X()+Tb*TDetectorReference.Y()+Tc*TDetectorReference.Z());
  //
  
  //
  Rotate3DX(angle);
  //
  
  return;
}

//____________________________________________________
void UNISLampWedgeDetector::RotateZ(Double_t angle)
{
  //
  //Rotation of the telescope reference frame to the input position
  TDetectorReference.RotateZ(angle);
  TDetectorTopReference.RotateZ(angle);
  TDetectorRightReference.RotateZ(angle);
  //
  
  //
  //calculation of the X Y versors
  TXversor=(TDetectorTopReference-TDetectorReference);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TDetectorRightReference-TDetectorReference);
  TYversor*=(1./TYversor.Mag());    
  //
  
  //
  //Calculation of the reference origin
  TDetectorNominalReference=TDetectorReference-TNominalDistanceBeamLine*TXversor;
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateZ(angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TDetectorReference.X()+Tb*TDetectorReference.Y()+Tc*TDetectorReference.Z());
  //
  
  //
  Rotate3DZ(angle);
  //
  
  return;
}

void UNISLampWedgeDetector::TranslateLongitudinal(Double_t z)
{
  //
  TDetectorReference+=z*TZlabversor;
  TDetectorTopReference+=z*TZlabversor;
  TDetectorRightReference+=z*TZlabversor;
  //
  
  //
  TranslateLongitudinal3D(z);
  //
  
  return;
}

void UNISLampWedgeDetector::Generate3D(double tilt, double phi_pos)
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
  fFrame       = new TGeoVolume("frame_volume",new TGeoTubeSeg(TNominalDistanceBeamLine, TNominalDistanceTopBeamLine, 0.05, -90-TFrameCoverageAngle*TMath::RadToDeg() , -90+TFrameCoverageAngle*TMath::RadToDeg()));
  fFrame->SetLineColor(kYellow+2);
  fDetector->AddNode(fFrame,0,new TGeoTranslation(0,TNominalDistanceBeamLine-TBottomFrame_distance,0));
  //
  
  //
  fStrip       = new TGeoVolume *[TStrips_number];
  //
  
  //
  for(Int_t i=0; i<TStrips_number; i++)
  {
    //
    fStrip[i] = new TGeoVolume(Form("strip_%02d_volume", i), new TGeoTubeSeg(TNominalDistanceBeamLine+TBottomFrame+TInter_width/2.+i*TStripTrue_width, TNominalDistanceBeamLine+TBottomFrame+TStripEffective_width+TInter_width/2.+i*TStripTrue_width, 0.1, 
                                                                             -90-TStripCoverageAngle[i]*TMath::RadToDeg() , -90+TStripCoverageAngle[i]*TMath::RadToDeg())); 
    fStrip[i]->SetLineColor(kGray);
    //
    fDetector->AddNode(fStrip[i],0,new TGeoTranslation(0,TNominalDistanceBeamLine-TBottomFrame_distance,0));
  }
  //
  
  //
  //We operate now a series of rotations to reach the final position
  //First: Rotation about the X-axis of a quantity (tilt)
  Rotate3DX(tilt);
  //Third: Rotation about the Z-axis of a quantity (phi)
  Rotate3DZ(phi_pos);
  //
}

void UNISLampWedgeDetector::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(x_angle*TMath::RadToDeg());
}

void UNISLampWedgeDetector::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISLampWedgeDetector::TranslateLongitudinal3D(Double_t z)
{
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(0.,0.,z));
}

// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
//____________________________________________________
Int_t UNISLampWedgeDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  Double_t t; /*t parameter for the parametric equation of the particle direction*/
  Double_t l,m,n; /*particle direction parameters*/
  l=TMath::Sin(theta_inc)*TMath::Cos(phi_inc);
  m=TMath::Sin(theta_inc)*TMath::Sin(phi_inc);
  n=TMath::Cos(theta_inc);

  //
  //solution of the direction/plane interception
  if((t=-(Ta*x0+Tb*y0+Tc*z0+Td)/(Ta*l+Tb*m+Tc*n))<0) return 0;
  //
  
  //
  //setting the true impact point
  TTrueImpactPoint.SetXYZ(x0+t*l,y0+t*m,z0+t*n);
  //  
  
  //
  //setting he impact point in the telescope frame
  TDetectorImpactPoint=TTrueImpactPoint-TDetectorNominalReference;
  double angle = TDetectorImpactPoint.Angle(TXversor);
  double distance = TDetectorImpactPoint.Mag();

  if (fabs(angle)<=TStripNominalCoverageAngle && distance>=TStripRadius[0]-TStripTrue_semi && distance<=TStripRadius[TStrips_number-1]+TStripTrue_semi) return 1;
  
  return 0;
}

// 
// Returns an absolute number identifying the pixel fired (in this case just strip)
// If the particle is not inside the active area -> return value = -1
//____________________________________________________
Int_t UNISLampWedgeDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  //Return -1 -> particle is not inside the detector
  
  double angle = TDetectorImpactPoint.Angle(TXversor);
  double distance = TDetectorImpactPoint.Mag();
  
  for(int i=0; i<TStrips_number; i++) {
    if(fabs(angle)<=TStripCoverageAngle[i] && fabs(distance-TStripRadius[i])<=TStripEffective_semi) return i;
  }
  
  return -1;
}

//draws the telescope on the X-Y plane
void UNISLampWedgeDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{}

//3D drawing function
void UNISLampWedgeDetector::Draw3D(Option_t * draw_opt) const
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

TVector3 UNISLampWedgeDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TTrueImpactPoint;
}

#ifdef GRAPHICAL_DEBUG
void UNISLampWedgeDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{}
#endif
