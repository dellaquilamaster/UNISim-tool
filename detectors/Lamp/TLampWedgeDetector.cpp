#include "TLampWedgeDetector.h"


//____________________________________________________
TLampWedgeDetector::TLampWedgeDetector(Double_t distance, Double_t phi_pos, Double_t tilt, Double_t bottom_frame_distance, Int_t N_Strips, Double_t strip_width, Double_t inter_width, Double_t nominal_radius, 
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
  TDetectorTopReference=TDetectorReference+TXlabversor;
  TDetectorRightReference=TDetectorReference+TYlabversor;
  //
  
  //
  //Rotation of the detector to the input tilt angle
  //
  RotateY(TTiltAngle);
  //

  //
  //Translation and rotation of the detector plane
  //
  TDetectorReference+=TDistanceBeamAxis*TZlabversor;
  TDetectorTopReference+=TDistanceBeamAxis*TZlabversor;
  TDetectorRightReference+=TDistanceBeamAxis*TZlabversor;
  TDetectorReference+=TBottomFrame_distance*TXlabversor;
  TDetectorTopReference+=TBottomFrame_distance*TXlabversor;
  TDetectorRightReference+=TBottomFrame_distance*TXlabversor;
  TDetectorReference.RotateZ(TAzimuthalAngle);
  TDetectorTopReference.RotateZ(TAzimuthalAngle);
  TDetectorRightReference.RotateZ(TAzimuthalAngle);
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
  OriginalDirection.RotateY(TTiltAngle);
  OriginalDirection.RotateZ(TAzimuthalAngle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TDetectorReference.X()+Tb*TDetectorReference.Y()+Tc*TDetectorReference.Z());
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
}

//____________________________________________________
TLampWedgeDetector::~TLampWedgeDetector()
{}

//____________________________________________________
void TLampWedgeDetector::RotateX(Double_t angle)
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
  
  return;
}

//____________________________________________________
void TLampWedgeDetector::RotateY(Double_t angle)
{
  //
  //Rotation of the telescope reference frame to the input position
  TDetectorReference.RotateY(angle);
  TDetectorTopReference.RotateY(angle);
  TDetectorRightReference.RotateY(angle);
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
  OriginalDirection.RotateY(angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TDetectorReference.X()+Tb*TDetectorReference.Y()+Tc*TDetectorReference.Z());
  //
  
  return;
}

// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
//____________________________________________________
Int_t TLampWedgeDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Int_t TLampWedgeDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
void TLampWedgeDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{}

//3D drawing function
void TLampWedgeDetector::Draw3D(Option_t * draw_opt) const
{
  //
  if(strstr(draw_opt,"SAME")==0 && strstr(draw_opt,"same")==0) {
    //
    //Creating the Graphics Manager
    TEveManager::Create();
    //
    
    //NOTE:
    //the Y-axis (first component of the vectors, horizontal axis) is considered with a "-" sign
    //for the rotation use the following scheme:
    // 3, 1 = rotation around X
    // 3, 2 = rotation around Y
    // 1, 2 = rotation around Z
    // Use RotatePF to rotate in the parent frame and RotateLF to rotate in the local frame
    //
    
    //
    //Generating Axes
    TEveGeoShape * TheAxes = new TEveGeoShape("TheAxes");
    gEve->AddElement(TheAxes);
    //
    TEveGeoShape * XaxisLine = new TEveGeoShape("XaxisLine");
    XaxisLine->SetShape(new TGeoTube(0, 0.2, 10));
    XaxisLine->SetMainColor(kGreen+1);
    XaxisLine->RefMainTrans().RotatePF(3, 2, 90*TMath::DegToRad());
    XaxisLine->RefMainTrans().Move3PF(0, 10, 0);
    XaxisLine->SetMainTransparency(50);
    TheAxes->AddElement(XaxisLine);
    TEveGeoShape * YaxisLine = new TEveGeoShape("YaxisLine");
    YaxisLine->SetShape(new TGeoTube(0, 0.2, 10));
    YaxisLine->SetMainColor(kRed);
    YaxisLine->RefMainTrans().RotatePF(3, 1, -90*TMath::DegToRad());
    YaxisLine->RefMainTrans().Move3PF(-10, 0, 0);
    YaxisLine->SetMainTransparency(50);
    TheAxes->AddElement(YaxisLine);
    TEveGeoShape * ZaxisLine = new TEveGeoShape("ZaxisLine");
    ZaxisLine->SetShape(new TGeoTube(0, 0.2, 10));
    ZaxisLine->SetMainColor(kBlue);
    ZaxisLine->RefMainTrans().Move3PF(0, 0, 10);
    ZaxisLine->SetMainTransparency(50);
    TheAxes->AddElement(ZaxisLine);
    TEveGeoShape * XaxisArrow = new TEveGeoShape("XaxisArrow");
    XaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    XaxisArrow->SetMainColor(kGreen+1);
    XaxisArrow->RefMainTrans().RotatePF(3, 2, 90*TMath::DegToRad());
    XaxisArrow->RefMainTrans().Move3PF(0, 20, 0);
    XaxisArrow->SetMainTransparency(50);
    TheAxes->AddElement(XaxisArrow);
    TEveGeoShape * YaxisArrow = new TEveGeoShape("YaxisArrow");
    YaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    YaxisArrow->SetMainColor(kRed);
    YaxisArrow->RefMainTrans().RotatePF(3, 1, -90*TMath::DegToRad());
    YaxisArrow->RefMainTrans().Move3PF(-20, 0, 0);
    YaxisArrow->SetMainTransparency(50);
    TheAxes->AddElement(YaxisArrow);
    TEveGeoShape * ZaxisArrow = new TEveGeoShape("ZaxisArrow");
    ZaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    ZaxisArrow->SetMainColor(kBlue);
    ZaxisArrow->RefMainTrans().Move3PF(0, 0, 20);
    ZaxisArrow->SetMainTransparency(50);
    TheAxes->AddElement(ZaxisArrow);
    TEveText * XaxisText = new TEveText("X");
    XaxisText->RefMainTrans().Move3PF(0, 21, 0);
    XaxisText->SetMainColor(kGreen+1);
    XaxisText->SetFontSize(20);
    XaxisText->SetLighting(kTRUE);
    TheAxes->AddElement(XaxisText);
    TEveText * YaxisText = new TEveText("Y");
    YaxisText->RefMainTrans().Move3PF(-21, 1, 0);
    YaxisText->SetMainColor(kRed);
    YaxisText->SetFontSize(20);
    YaxisText->SetLighting(kTRUE);
    TheAxes->AddElement(YaxisText);
    TEveText * ZaxisText = new TEveText("Z");
    ZaxisText->RefMainTrans().Move3PF(0, 1, 21);
    ZaxisText->SetMainColor(kBlue);
    ZaxisText->SetFontSize(20);
    ZaxisText->SetLighting(kTRUE);
    TheAxes->AddElement(ZaxisText);
    //
  }
  //

  //
  //Generating detector frame
  TEveGeoShape * DetectorFrame = new TEveGeoShape("DetectorFrame");
  //TGeoTubeSeg: ----------------------------> name -----------> rmin --------------> rmax ------------------------- > dz ---> phi1 --------------------------------------> phi2
  TGeoShape * TheFrameShape = new TGeoTubeSeg("TheFrameShape", TNominalDistanceBeamLine, TNominalDistanceTopBeamLine, 0.05, -TFrameCoverageAngle*TMath::RadToDeg() , TFrameCoverageAngle*TMath::RadToDeg());
  DetectorFrame->SetShape(TheFrameShape);
  DetectorFrame->SetMainColor(kYellow+2);
  
  /*
  // TGeoTrap: -----------------------------> name ------------> height -------------------------->  lower base ---------------> thickness ------> upper base -----------------------> thickness
  TGeoShape * UpperFrameShape = new TGeoTrap ("UpperFrameShape", (THeight-THeightPrime)/2., 0., 0., TTopBase/2.+TLateralFrame/2., 0.05, 0.05, 0., TDetectorSemiWidth+TLateralFrame/2., 0.05, 0.05, 0.); //Creating shape for the upper wedge
  TGeoShape * LowerFrameShape = new TGeoTrap ("LowerFrameShape", THeightPrime/2., 0., 0., TBottomBase/2.+TLateralFrame/2., 0.05, 0.05, 0., TDetectorSemiWidth+TLateralFrame/2., 0.05, 0.05, 0.); //Creating shape for the lower wedge
  TGeoRotation *RotUpper = new TGeoRotation("RotUpper",180,0.,90.,0.,90,270); //(Identity matrix is 90, 0, 90, 90, 0, 0)
  TGeoRotation *RotLower = new TGeoRotation("RotLower",180,0.,90.,180.,90,90); //(Identity matrix is 90, 0, 90, 90, 0, 0)
  TGeoCombiTrans * MatrixUpper = new TGeoCombiTrans("MatrixUpper",0.,(THeight-THeightPrime)/2.+THeightPrime,0,RotUpper); //Transformation matrix for the upper wedge
  TGeoCombiTrans * MatrixLower = new TGeoCombiTrans("MatrixLower",0.,THeightPrime/2.,0,RotLower); //Transformation matrix for the lower wedge
  MatrixUpper->RegisterYourself(); //Registering Transformation matrix of the upper wedge
  MatrixLower->RegisterYourself(); //Registering Transformation matrix of the lower wedge
  DetectorFrame->SetShape(new TGeoCompositeShape("TheFrame", "(UpperFrameShape:MatrixUpper+LowerFrameShape:MatrixLower)"));
  //
  */
    
  //
  //Rotation of the frame to the input position
  DetectorFrame->RefMainTrans().RotatePF(1, 2, TMath::Pi()/2.);
  DetectorFrame->RefMainTrans().Move3PF(0., -TNominalDistanceBeamLine,0.);
  DetectorFrame->RefMainTrans().RotatePF(3, 2, TTiltAngle);
  DetectorFrame->RefMainTrans().Move3PF(0., TBottomFrame_distance, TDistanceBeamAxis);
  DetectorFrame->RefMainTrans().RotatePF(1, 2, TAzimuthalAngle);
  //
  
  //
  //Adding detector frame to TEveManager
  gEve->AddElement(DetectorFrame);
  //
  
  //
  //Generating strips
  TEveGeoShape * Strip[TStrips_number];
  for(int i=0; i<TStrips_number; i++)
  {
    //TGeoTubeSeg ------------->
    Strip[i] = new TEveGeoShape(Form("Strip_%02d",i));
    Strip[i]->SetShape(new TGeoTubeSeg(Form("Strip_%02d_Shape", i), TNominalDistanceBeamLine+TBottomFrame+TInter_width/2.+i*TStripTrue_width, TNominalDistanceBeamLine+TBottomFrame+TStripEffective_width+TInter_width/2.+i*TStripTrue_width, 0.1, -TStripCoverageAngle[i]*TMath::RadToDeg() , TStripCoverageAngle[i]*TMath::RadToDeg()));
    Strip[i]->SetMainColor(kGray);
    Strip[i]->RefMainTrans().RotatePF(1, 2, TMath::Pi()/2.);
    Strip[i]->RefMainTrans().Move3PF(0., -TNominalDistanceBeamLine,0.);
    Strip[i]->RefMainTrans().RotatePF(3, 2, TTiltAngle);
    Strip[i]->RefMainTrans().Move3PF(0., TBottomFrame_distance, TDistanceBeamAxis);
    Strip[i]->RefMainTrans().RotatePF(1, 2, TAzimuthalAngle);
    DetectorFrame->AddElement(Strip[i]);
  }
  //
  
  // Drawing3D  
  gEve->Redraw3D(kTRUE);
  //
}

TVector3 TLampWedgeDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TTrueImpactPoint;
}

#ifdef GRAPHICAL_DEBUG
void TLampWedgeDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{}
#endif
