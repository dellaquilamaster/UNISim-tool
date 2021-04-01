#include "UNISLampWedgeMMMDetector.h"


//____________________________________________________
UNISLampWedgeMMMDetector::UNISLampWedgeMMMDetector(Double_t distance, Double_t phi_pos, Double_t tilt, Double_t bottom_frame_distance, Double_t inter_width) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TDistanceBeamAxis(distance),
TAzimuthalAngle(phi_pos),
TTiltAngle(tilt),
TBottomFrame_distance(bottom_frame_distance),
TAnnularStrips_number(16),                                                  
TRadialStrips_number(8),                                                  
TAnnularStripTrue_width(0.6410),                                             
TAnnularStripTrue_semi(0.5*TAnnularStripTrue_width),                                                                        
TInter_width(inter_width),   
TAnnularStripEffective_width(TAnnularStripTrue_width-TInter_width),
TAnnularStripEffective_semi(TAnnularStripEffective_width/2.),
TInnerNominal_radius(6.520),
TStripNominalCoverageAngle(56./2.*TMath::DegToRad()), //TO BE FIXED
TAnnularStripCoverageAngle(new Double_t[TAnnularStrips_number]),
TRadialStripMinimumEffectiveAngle(new Double_t[TRadialStrips_number]),
TRadialStripMaximumEffectiveAngle(new Double_t[TRadialStrips_number]),
TFrameCoverageAngle(60./2.*TMath::DegToRad()),                                                    
TBottomFrame(1.),
TTopFrame_distance(TBottomFrame_distance+TBottomFrame+TAnnularStripTrue_width*TAnnularStrips_number+TBottomFrame),
TNominalDistanceBeamLine(TInnerNominal_radius-TBottomFrame),
TNominalDistanceTopBeamLine(TNominalDistanceBeamLine+2*TBottomFrame+TAnnularStrips_number*TAnnularStripTrue_width),
TAnnularStripRadius(new Double_t[TAnnularStrips_number])
{
  //
  //Placing the detector plane orthogonal to the beam axis and with its bottom center touching the beam line
  //The plane is identified by the three reference points (origin on the beam-line and 2 additional points)
  TDetectorReference.SetXYZ(0.,0.,0.);
  TDetectorTopReference=TDetectorReference+TYlabversor;
  TDetectorRightReference=TDetectorReference-TXlabversor;
  //
  
  //
  //Calculation of the reference origin
  TDetectorNominalReference=TDetectorReference-TNominalDistanceBeamLine*TXversor;
  //
  
  //
  //Generating strips
  for(int i=0; i<TAnnularStrips_number; i++) {
    //
    TAnnularStripRadius[i]=TInnerNominal_radius+(2*i+1)*TAnnularStripTrue_semi;
    TAnnularStripCoverageAngle[i]=TStripNominalCoverageAngle;
    //
  }
  for(int i=0; i<TRadialStrips_number; i++) {
    //
    TRadialStripMinimumEffectiveAngle[i]=-TStripNominalCoverageAngle+i*TStripNominalCoverageAngle*2/TRadialStrips_number+asin(TInter_width/(2.*TInnerNominal_radius));
    TRadialStripMaximumEffectiveAngle[i]=-TStripNominalCoverageAngle+(i+1)*TStripNominalCoverageAngle*2./TRadialStrips_number-asin(TInter_width/(2.*TInnerNominal_radius));
    //
  }
  //
  
  
  //
  //calculation of the X Y versors identifying the detector plane
  TXversor=(TDetectorTopReference-TDetectorReference);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TDetectorRightReference-TDetectorReference);
  TYversor*=(1./TYversor.Mag());
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
  //Translate the detector to the desired distance from the beamline
  Translate(TBottomFrame_distance, TYlabversor);
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
UNISLampWedgeMMMDetector::~UNISLampWedgeMMMDetector()
{}

//____________________________________________________
void UNISLampWedgeMMMDetector::RotateX(Double_t angle)
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
void UNISLampWedgeMMMDetector::RotateZ(Double_t angle)
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


void UNISLampWedgeMMMDetector::TranslateLongitudinal(Double_t z)
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

void UNISLampWedgeMMMDetector::Translate(Double_t quantity, TVector3 vers)
{
  //
  TVector3 TheVersor(vers.Unit());
  //
  
  //
  //Moving reference point according to the bottom frame displacement
  TDetectorReference+=TBottomFrame_distance*TheVersor;
  TDetectorTopReference+=TBottomFrame_distance*TheVersor; 
  TDetectorRightReference+=TBottomFrame_distance*TheVersor;
  //
  
  //
  //Calculation of the reference origin
  TDetectorNominalReference+=quantity*TheVersor;
  //
  
  //
  Translate3D(quantity, TheVersor);
  //
  
  return;
}

void UNISLampWedgeMMMDetector::Generate3D(double tilt, double phi_pos)
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
  fFrame->SetLineColor(kGreen+2);
  fDetector->AddNode(fFrame,0,new TGeoTranslation(0,TNominalDistanceBeamLine,0));  
  //
  
  //
  fStripAnnular= new TGeoVolume *[TAnnularStrips_number];
  fStripRadial = new TGeoVolume *[TRadialStrips_number];
  //
  
  //
  for(int i=0; i<TAnnularStrips_number; i++)
  {
    fStripAnnular[i] = new TGeoVolume(Form("strip_annular_%02d_volume", i), new TGeoTubeSeg(TNominalDistanceBeamLine+TBottomFrame+TInter_width/2.+i*TAnnularStripTrue_width, TNominalDistanceBeamLine+TBottomFrame+TAnnularStripEffective_width+TInter_width/2.+i*TAnnularStripTrue_width, 0.1, -90-TAnnularStripCoverageAngle[i]*TMath::RadToDeg() , -90+TAnnularStripCoverageAngle[i]*TMath::RadToDeg())); 
    //
    fStripAnnular[i]->SetLineColor(kGray);
    //
    fDetector->AddNode(fStripAnnular[i],0,new TGeoTranslation(0,TNominalDistanceBeamLine,-0.1));
    //
  }
  for(int i=0; i<TRadialStrips_number; i++)
  {
    //
    fStripRadial[i] = new TGeoVolume(Form("strip_radial_%02d_volume", i), new TGeoTubeSeg(TNominalDistanceBeamLine+TBottomFrame+TInter_width/2., TNominalDistanceBeamLine+TBottomFrame+TAnnularStripTrue_width*TAnnularStrips_number-TInter_width/2., 0.1, -90+TRadialStripMinimumEffectiveAngle[i]*TMath::RadToDeg() , -90+TRadialStripMaximumEffectiveAngle[i]*TMath::RadToDeg()));
    //
    fStripRadial[i]->SetLineColor(kGray);
    //
    fDetector->AddNode(fStripRadial[i],0,new TGeoTranslation(0,TNominalDistanceBeamLine,0.1));
    //
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

void UNISLampWedgeMMMDetector::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(-x_angle*TMath::RadToDeg());
  //NOTE: x-axis is swaped compared to that used in the calculations, so I need to use -x_angle 
}

void UNISLampWedgeMMMDetector::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISLampWedgeMMMDetector::TranslateLongitudinal3D(Double_t z)
{
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(0.,0.,z));
}

void UNISLampWedgeMMMDetector::Translate3D(Double_t quantity, TVector3 vers)
{
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(quantity*vers.X(),-quantity*vers.Y(),quantity*vers.Z()));
  //NOTE: y-axis is swaped compared to the axis used in the calculations
}

// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
//____________________________________________________
Int_t UNISLampWedgeMMMDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
    
  if (fabs(angle)<=TStripNominalCoverageAngle && distance>=TAnnularStripRadius[0]-TAnnularStripTrue_semi && distance<=TAnnularStripRadius[TAnnularStrips_number-1]+TAnnularStripTrue_semi) return 1;
  
  return 0;
}

// 
// Returns an absolute number identifying the pixel fired (in this case just strip)
// If the particle is not inside the active area -> return value = -1
//____________________________________________________
Int_t UNISLampWedgeMMMDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  //Return -1 -> particle is not inside the detector
  
  double angle = TDetectorImpactPoint.Angle(TXversor)*TDetectorImpactPoint.Dot(TYversor)/std::fabs(TDetectorImpactPoint.Dot(TYversor));
  double distance = TDetectorImpactPoint.Mag();
  
  int annular_strip=-1;
  int radial_strip=-1;
  
//   printf("%f %f %f\n", TDetectorImpactPoint.X(), TDetectorImpactPoint.Y(), TDetectorImpactPoint.Z());
//   printf("angle=%f distance=%f\n", angle*TMath::RadToDeg(),distance);
  
  //
  for(int i=0; i<TAnnularStrips_number; i++) {
    if(fabs(angle)<=TAnnularStripCoverageAngle[i] && fabs(distance-TAnnularStripRadius[i])<=TAnnularStripEffective_semi) annular_strip=i;
  }
  for(int i=0; i<TRadialStrips_number; i++) {
    if(angle>=TRadialStripMinimumEffectiveAngle[i] && angle<=TRadialStripMaximumEffectiveAngle[i]) radial_strip=i;
  }
  //
  
//   printf("%d %d -> pixel=%d\n", annular_strip, radial_strip, annular_strip*TRadialStrips_number+radial_strip);
  
  //
  if(annular_strip>=0 && radial_strip>=0) {
    return annular_strip*TRadialStrips_number+radial_strip;
  }
  //
  
  return -1;
}

//draws the telescope on the X-Y plane
void UNISLampWedgeMMMDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{}

//3D drawing function
void UNISLampWedgeMMMDetector::Draw3D(Option_t * draw_opt) const
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

TVector3 UNISLampWedgeMMMDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TTrueImpactPoint;
}

#ifdef GRAPHICAL_DEBUG
void UNISLampWedgeMMMDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{}
#endif
