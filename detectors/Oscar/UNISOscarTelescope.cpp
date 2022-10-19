#include "UNISOscarTelescope.h"

//____________________________________________________
UNISOscarTelescope::UNISOscarTelescope(Double_t distance, Double_t theta_pos, Double_t phi_pos, Option_t * opt) :
fNumPads(16),
fNumStrips(16),
fPadRowsColumns(sqrt(fNumPads)),
fPadWidth(1.),
fCollimatorWidth(0.8),
fPadSemi(fPadWidth/2.),
fCollimatorSemi(fCollimatorWidth/2.),
fPadFrameWidth(0.14),
fPadBottomContactsWidth(0.18),
fPhotoDiodeWidth(2*fPadFrameWidth+fPadWidth),
fPhotoDiodeHeight(2*fPadFrameWidth+fPadWidth+fPadBottomContactsWidth),
fFrameWidth((fPhotoDiodeWidth*fPadRowsColumns)+1),
fFrameHeight((fPhotoDiodeHeight*fPadRowsColumns)+2.),
fIsStrip(true),
fIsCollimator(false),
fStripPadDistance(0.5),
fXlabversor(1,0,0),
fYlabversor(0,1,0),
fZlabversor(0,0,1),
fCenter(0,0,0),
fLabImpactPoint(0.,0.,0.),
fFrameImpactPoint(0.,0.,0.),
fPads(new UNISSiliconPhotoDiode * [fNumPads])
{
  //
  if(std::string(opt).find("pads")!=std::string::npos) fIsStrip=false;
  if(std::string(opt).find("col")!=std::string::npos) fIsCollimator=true;
  
  //
  //Telescope's corners (telescope at zero-distance)
  fTopLeftCorner.SetXYZ(fFrameWidth/2.,fFrameHeight/2.,0.);
  fTopRightCorner.SetXYZ(-fFrameWidth/2.,fFrameHeight/2.,0.);
  fBottomLeftCorner.SetXYZ(fFrameWidth/2.,-fFrameHeight/2.,0.);
  fBottomRightCorner.SetXYZ(-fFrameWidth/2.,-fFrameHeight/2.,0.);
  //
  
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftCorner-fBottomLeftCorner);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightCorner-fTopLeftCorner);
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
  //Creation of a telescope along the beam axis at zero-distance from the target
  //
  //pads
  for(int row=0; row<fPadRowsColumns; row++) {
    for(int col=0; col<fPadRowsColumns; col++) {
      const double y_coordinate = (2*std::fabs((fPadRowsColumns-1)/2.-row)*(fPhotoDiodeWidth/2.)+(std::fabs((fPadRowsColumns-1)/2.-row)-0.5)*(fPadBottomContactsWidth));
      fPads[row*fPadRowsColumns+col] = new UNISSiliconPhotoDiode((fPadRowsColumns-1)*fPhotoDiodeWidth/2.-col*fPhotoDiodeWidth, //x-coordinate (horizontal towards left, as seen from the beam)
                                                                 (row < fPadRowsColumns/2 ? y_coordinate : -y_coordinate), //y-coordinate (vertical towards the top, as seen from the beam)
                                                                 0., //z-coordinate
                                                                 0,0,row<fPadRowsColumns/2 ? TMath::Pi() : 0, //tilt-angles //inserting pad with only z-tilt
                                                                 (fIsCollimator ? fCollimatorWidth : fPadWidth)); //effective width of the collimator
    }
  }
  //
  //strips
  if(fIsStrip) {
    fStrip = new UNISStripSingleSidedDetector(0.,0.,-0.3,0.,0.,16,0.312,0.01,0.4,0.);
  }
  //
  
  //
  Generate3D();
  //
  
  
  //
  Translate(0,0,distance);
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

//____________________________________________________
UNISOscarTelescope::UNISOscarTelescope(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y, Option_t * opt) :
UNISOscarTelescope(0,0,0,opt)
{
  //NOTE: to be implemented
}

//____________________________________________________
UNISOscarTelescope::~UNISOscarTelescope()
{}

//____________________________________________________
void UNISOscarTelescope::RotateX(Double_t angle)
{
  //
  fCenter.RotateX(angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftCorner   .RotateX(angle); 
  fTopRightCorner  .RotateX(angle); 
  fBottomLeftCorner.RotateX(angle);
  fBottomRightCorner.RotateX(angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftCorner-fBottomLeftCorner);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightCorner-fTopLeftCorner);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateX(angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  //Rotating pads
  for(int pad=0; pad<fNumPads; pad++) {
    fPads[pad]->RotateX(angle); 
  }
  //
  
  //
  //Rotating strips
  if(fIsStrip) {
    fStrip->RotateX(angle); 
  }
  //
  
  //
  //Rotating 3D frame
  fDetectorMatrix->RotateX(angle*TMath::RadToDeg());
  //
}

//____________________________________________________
void UNISOscarTelescope::RotateY(Double_t angle)
{
  //
  fCenter.RotateY(angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftCorner   .RotateY(angle); 
  fTopRightCorner  .RotateY(angle); 
  fBottomLeftCorner.RotateY(angle);
  fBottomRightCorner.RotateY(angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftCorner-fBottomLeftCorner);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightCorner-fTopLeftCorner);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateY(angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  //Rotating pads
  for(int pad=0; pad<fNumPads; pad++) {
    fPads[pad]->RotateY(angle); 
  }
  //
  
  //
  //Rotating strips
  if(fIsStrip) {
    fStrip->RotateY(angle); 
  }
  //
  
  //
  //Rotating 3D frame
  fDetectorMatrix->RotateY(angle*TMath::RadToDeg());
  //
}

//____________________________________________________
void UNISOscarTelescope::RotateZ(Double_t angle)
{
  //
  fCenter.RotateZ(angle);
  //
  
  //
  //rotation of the corners to the input position
  fTopLeftCorner   .RotateZ(angle); 
  fTopRightCorner  .RotateZ(angle); 
  fBottomLeftCorner.RotateZ(angle);
  fBottomRightCorner.RotateZ(angle);
  //
  //calculation of the X Y versors
  fYversor=(fTopLeftCorner-fBottomLeftCorner);
  fYversor*=(1./fYversor.Mag());
  fXversor=(fTopRightCorner-fTopLeftCorner);
  fXversor*=(1./fXversor.Mag());  
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(fa,fb,fc);
  OriginalDirection.RotateZ(angle);
  fa=OriginalDirection.X()/OriginalDirection.Mag();
  fb=OriginalDirection.Y()/OriginalDirection.Mag();
  fc=OriginalDirection.Z()/OriginalDirection.Mag();
  fd=-(fa*fCenter.X()+fb*fCenter.Y()+fc*fCenter.Z());
  //
  
  //
  //Rotating pads
  for(int pad=0; pad<fNumPads; pad++) {
    fPads[pad]->RotateZ(angle); 
  }
  //
  
  //
  //Rotating strips
  if(fIsStrip) {
    fStrip->RotateZ(angle); 
  }
  //
  
  //
  //Rotating 3D frame
  fDetectorMatrix->RotateZ(angle*TMath::RadToDeg());
  //
}

//____________________________________________________
void UNISOscarTelescope::Translate(Double_t x, Double_t y, Double_t z)
{
  //
  fCenter+=TVector3(x,y,z);
  //  
  
  //
  //Translating pads
  for(int pad=0; pad<fNumPads; pad++) {
    fPads[pad]->Translate(x,y,z); 
  }
  //
  
  //
  //Translating strips
  if(fIsStrip) {
    fStrip->Translate(x,y,z);
  }
  //
  
  //
  //Translating 3D frame
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(x,y,z));
  //
}

//
// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
Int_t UNISOscarTelescope::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
    
  // setting the impact point in the telescope frame
  fFrameImpactPoint=fLabImpactPoint-fCenter;
  const double ImpactXprime = fFrameImpactPoint.Dot(fXversor);
  const double ImpactYprime = fFrameImpactPoint.Dot(fYversor); 
  
  //
  //Inside condition
  if(std::fabs(ImpactXprime)<=fFrameWidth/2. && std::fabs(ImpactYprime)<=fFrameHeight/2.)
  {
    return 1;
  }
  //
  return 0;
}

// 
// Returns 0 if the particle is inside the active area:
// If the particle is not inside the active area -> return value = -1
Int_t UNISOscarTelescope::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  //
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  /*If not the particle is inside the telescope*/
  //
    
  int PixelPad=-1;
  int PixelStrip;
    
  //
  //Check if the particle is inside the active area of a pad
  for(int pad=0; pad<fNumPads; pad++) {
    if(fPads[pad]->GetPixel(theta_inc, phi_inc, x0, y0, z0)!=-1) {
      PixelPad=pad;
      break;
    }
  }
  //
  
  //
  if(PixelPad==-1) return -1; //particle not inside a pad
  //
  
  //
  //Check if the particle is inside the active area of a strip (if strip is present)
  if(fIsStrip) {
    PixelStrip=fStrip->GetPixel(theta_inc, phi_inc, x0, y0, z0);
    if(PixelStrip==-1) return -1; //particle not inside a strip (and strip is installed)
  }
  //
  
  //
  //Calculation of pixel
  int Pixel=fIsStrip ? (PixelPad)*fPadRowsColumns+(PixelStrip%fPadRowsColumns) : PixelPad;
  //
  
  //particle not inside the effective area
  return Pixel;
}

//____________________________________________________
void UNISOscarTelescope::Generate3D()
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
  fDetector = new TGeoVolumeAssembly("DetectorFrame");
  fDetectorMatrix = new TGeoHMatrix("DetectorTransformationMatrix");
  //
   
  //
  fFramePads      = new TGeoVolume("frame_pads_volume",new TGeoBBox(fFrameWidth/2., fFrameHeight/2., 0.08));
  fFramePreAmps   = new TGeoVolume("frame_preamps_volume",new TGeoBBox((fPhotoDiodeWidth*fPadRowsColumns)/2.+0.5, (fPhotoDiodeHeight*fPadRowsColumns)/2.+1., 0.08));
  fPreAmps        = new TGeoVolume("preamps_volume",new TGeoBBox(1., 0.15, 1.));
  //
  TGeoBBox * collimator_frame = new TGeoBBox("collimator_frame", fFrameWidth/2., fFrameHeight/2., 0.04);
  TGeoBBox * collimator_hole = new TGeoBBox("collimator_hole", fCollimatorWidth/2., fCollimatorWidth/2., 0.041);
  TGeoBBox * collimator_bottom_hole = new TGeoBBox("collimator_bottom_hole", fFrameWidth/2.-0.8, 0.8, 0.041);
  TGeoTranslation ** collimator_hole_trans = new TGeoTranslation * [fNumPads];
  std::string compisite_shape_string("collimator_frame - ");
  for(int pad=0; pad<fNumPads; pad++) {
    const int row=pad/4;
    const int col=pad%4;
    const double y_coordinate = (2*std::fabs((fPadRowsColumns-1)/2.-row)*(fPhotoDiodeWidth/2.)+(std::fabs((fPadRowsColumns-1)/2.-row)-0.5)*(fPadBottomContactsWidth));
    collimator_hole_trans[pad]=new TGeoTranslation(Form("collimator_hole_trans_%02d", pad), -(fPadRowsColumns-1)*fPhotoDiodeWidth/2.+col*fPhotoDiodeWidth, (row < fPadRowsColumns/2 ? y_coordinate : -y_coordinate)+0.5,0.);
    collimator_hole_trans[pad]->SetName(Form("collimator_hole_trans_%02d",pad));
    collimator_hole_trans[pad]->RegisterYourself();
    compisite_shape_string.append(Form("collimator_hole:collimator_hole_trans_%02d %s ", pad, pad==fNumPads-1 ? "" : "-"));
  }
  TGeoTranslation * collimator_bottom_hole_trans = new TGeoTranslation(0., -fFrameHeight/2.+0.25, 0.);
  collimator_bottom_hole_trans->SetName("collimator_bottom_hole_trans");
  collimator_bottom_hole_trans->RegisterYourself();
  compisite_shape_string.append(" - collimator_bottom_hole:collimator_bottom_hole_trans");
  //
  fCollimator = new TGeoVolume("collimator", new TGeoCompositeShape("collimator_shape", compisite_shape_string.c_str()));
  fFramePads->SetLineColor(kGreen+2);
  fFramePreAmps->SetLineColor(kGreen+2);
  fPreAmps->SetLineColor(kMagenta);
  fCollimator->SetLineColor(kOrange+1);
  //
  
  //
  //Adding to mother volume
  fDetector->AddNode(fFramePads,0,new TGeoTranslation(0, -0.5, 0.)); //distance 0. cm
  fDetector->AddNode(fFramePreAmps,0,new TGeoTranslation(0, -0.5, 0.4)); //distance 0.4 cm
  if(fIsCollimator) fDetector->AddNode(fCollimator,0,new TGeoTranslation(0, -0.5, -0.38)); //distance -0.3 cm
  for(int pad=0; pad<fNumPads/2.; pad++) {
    fDetector->AddNode(fPreAmps,pad,new TGeoTranslation(fFrameWidth/4., fFrameHeight/2.*0.6-pad*0.5, 0.4+0.08+1)); //on the preamp board (0.4 cm)
    fDetector->AddNode(fPreAmps,pad,new TGeoTranslation(-fFrameWidth/4., fFrameHeight/2.*0.6-pad*0.5, 0.4+0.08+1)); //on the preamp board (0.4 cm)
  }
  //
  
  //
  fDetectorMatrix->MultiplyLeft(new TGeoTranslation(0,0,0.2));
  //
}

//____________________________________________________
void UNISOscarTelescope::Draw3D(Option_t * draw_opt) const
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
  
  //
  //Drawing pads
  for(int pad=0; pad<fNumPads; pad++) {
    fPads[pad]->Draw3D(Form("%s SAME", draw_opt)); 
  }
  //
  
  //
  //Translating strips
  if(fIsStrip) {
    fStrip->Draw3D(Form("%s SAME", draw_opt));
  }
  //
  //
}

//____________________________________________________
void UNISOscarTelescope::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return; 
}

//____________________________________________________
TVector3 UNISOscarTelescope::GetDetectorCenter()
{
  return fCenter; 
}

//____________________________________________________
TVector3 UNISOscarTelescope::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return fLabImpactPoint;
}
