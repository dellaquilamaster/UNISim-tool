#include "UNISFaziaQuartet.h"

// standard constructor
UNISFaziaQuartet::UNISFaziaQuartet(Double_t theta_pos, Double_t phi_pos, Double_t displacement, Double_t pad_width, Double_t frame_width) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TPads_number(4),
TRowColumn(2),
TNominalDistance(100),
TDisplacement(displacement),
TNominalTheta(theta_pos),
TNominalPhi(phi_pos),
TPadEffective_width(pad_width),
TPadEffective_semi(0.5*TPadEffective_width),
TFrame_width(frame_width),
TPadTrue_width(TPadEffective_width+2*TFrame_width),
TPadTrue_semi(0.5*TPadTrue_width),
TTelescopeTrue_semi(TPadTrue_semi*TRowColumn),
TTrueImpactPoint(0.,0.,0.),
TTelescopeImpactPoint(0.,0.,0.)
{
  //
  //Placing the detector at 0 degrees and distance equal to the nominal one (100 cm)
  //
  //Telescope corners
  TTopLeftCorner.SetXYZ(TTelescopeTrue_semi,-TTelescopeTrue_semi,TNominalDistance);
  TTopRightCorner.SetXYZ(TTelescopeTrue_semi,TTelescopeTrue_semi,TNominalDistance);
  TBottomLeftCorner.SetXYZ(-TTelescopeTrue_semi,-TTelescopeTrue_semi,TNominalDistance);
  TTopLeftXCorner= TTelescopeTrue_semi;
  TTopLeftYCorner= TTelescopeTrue_semi;
  /* vectors allocations*/
  TCenters       =new TVector3*[TRowColumn];
  TCentersXprime =new Double_t*[TRowColumn];
  TCentersYprime =new Double_t*[TRowColumn];
  for(Int_t i = 0; i<TRowColumn; i++)
  {
    TCenters[i]       = new TVector3[TRowColumn]; 
    TCentersXprime[i] = new Double_t[TRowColumn]; 
    TCentersYprime[i] = new Double_t[TRowColumn]; 
  }
  //
  //Telescope center
  TTelescopeCenter.SetXYZ(0.,0.,TNominalDistance);
  //
  //Telescope reference frame versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
  // placement of the entire cluster of pixels on a plane orthogonal to z-axis at the input TNominalDistance
  // here i represents the pad row and j represents the pad column
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].SetXYZ((TRowColumn-(2*j+1))*TPadTrue_semi,-(TRowColumn-(2*i+1))*TPadTrue_semi,TNominalDistance);
    }
  }
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
  Generate3D(0.,0.);
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
  
  //
  //Displacement of the cluster
  Translate(0,0,TDisplacement);
  //
}

void UNISFaziaQuartet::RotateZ(Double_t z_angle)
{
  //Rotation of the telescope center
  TTelescopeCenter.RotateZ(z_angle);
  
  //Rotation of the corners
  TTopLeftCorner   .RotateZ(z_angle); 
  TTopRightCorner  .RotateZ(z_angle); 
  TBottomLeftCorner.RotateZ(z_angle);
  
  //Rotation of pixels
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateZ(z_angle);
    }
  }
  
  //Calculation of detector reference frame vectors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag()); 
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateZ(z_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the pad, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
    }
  }
  //
  
  //
  Rotate3DZ(z_angle);
  //
  
  return;
}

void UNISFaziaQuartet::Translate(Double_t x, Double_t y, Double_t z)
{
  //
  TVector3 TranslationVector (x, y, z);
  
  //Translation of the telescope center
  TTelescopeCenter+=TranslationVector;
  
  //Translation of the corners
  TTopLeftCorner+=TranslationVector;
  TTopRightCorner+=TranslationVector;
  TBottomLeftCorner+=TranslationVector;
  
  //Translation of pixels
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j]+=TranslationVector;
    }
  }
  
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
  
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the pad, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
    }
  }
  //
  
  //
  Translate3D(x,y,z);
  //
  
  return;
}

void UNISFaziaQuartet::RotateX(Double_t x_angle)
{
 //Rotation of the telescope center
  TTelescopeCenter.RotateX(x_angle);
  
  //Rotation of the corners
  TTopLeftCorner   .RotateX(x_angle); 
  TTopRightCorner  .RotateX(x_angle); 
  TBottomLeftCorner.RotateX(x_angle);
  
  //Rotation of pixels
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateX(x_angle);
    }
  }
  
  //Calculation of detector reference frame vectors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag()); 
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateX(x_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the pad, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
    }
  }
  //
  
  //
  Rotate3DX(x_angle);
  //
  
  return;
}

void UNISFaziaQuartet::RotateY(Double_t y_angle)
{
  //Rotation of the telescope center
  TTelescopeCenter.RotateY(y_angle);
  
  //Rotation of the corners
  TTopLeftCorner   .RotateY(y_angle); 
  TTopRightCorner  .RotateY(y_angle); 
  TBottomLeftCorner.RotateY(y_angle);
  
  //Rotation of pixels
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateY(y_angle);
    }
  }
  
  //Calculation of detector reference frame vectors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag()); 
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateY(y_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the pad, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
    }
  }
  //
  
  //
  Rotate3DY(y_angle);
  //
  
  return;
}

void UNISFaziaQuartet::Generate3D(double theta_pos, double phi_pos)
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
  fPad         = new TGeoVolume("pad_volume",new TGeoBBox(TPadEffective_semi, TPadEffective_semi, 0.1));
  fTopFrame    = new TGeoVolume("top_frame_volume",new TGeoBBox(TPadTrue_semi,TFrame_width/2., 0.2));
  fBottomFrame = new TGeoVolume("bottom_frame_volume",new TGeoBBox(TPadTrue_semi,TFrame_width/2., 0.2));
  fLeftFrame   = new TGeoVolume("left_frame_volume",new TGeoBBox(TFrame_width/2.,TPadTrue_semi, 0.2, 0));
  fRightFrame  = new TGeoVolume("right_frame_volume",new TGeoBBox(TFrame_width/2.,TPadTrue_semi, 0.2, 0));
  fPad->SetLineColor(kGray);
  fTopFrame->SetLineColor(kYellow+2);
  fBottomFrame->SetLineColor(kYellow+2);
  fLeftFrame->SetLineColor(kYellow+2);
  fRightFrame->SetLineColor(kYellow+2);
  //
  
  //
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      fDetector->AddNode(fPad,i*TRowColumn+j,new TGeoTranslation(-(TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,0));
      fDetector->AddNode(fTopFrame,i*TRowColumn+j,new TGeoTranslation(-(TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,0));
      fDetector->AddNode(fBottomFrame,i*TRowColumn+j,new TGeoTranslation(-(TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,0));
      fDetector->AddNode(fLeftFrame,i*TRowColumn+j,new TGeoTranslation(-(TRowColumn-(2*i+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0));
      fDetector->AddNode(fRightFrame,i*TRowColumn+j,new TGeoTranslation(-(TRowColumn-(2*i+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0));
    }
  }
  //
  
  //
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(0,0,TNominalDistance));
  //
  
  //
  //Drawing decoratively a block of CsIs
  const double fCsICrystalLength=10.; // cm
  fCsICrystal = new TGeoVolume("crystal_volume",new TGeoPara(TPadEffective_semi, TPadEffective_semi, fCsICrystalLength/2., 0, 0, 0));
  fCsICrystal->SetLineColor(kBlue-4);
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      fDetector->AddNode(fCsICrystal,i*TRowColumn+j,new TGeoTranslation((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,fCsICrystalLength/2.+0.1));
    }
  }
  //
  
  //
  //We operate now a series of rotations to reach the final position
  //First: Rotation about the Z-axis of a quantity (-phi)
  Rotate3DZ(-phi_pos-180*TMath::DegToRad()); 
  //Second: Rotation about the X-axis of a quantity (theta)
  Rotate3DX(theta_pos);
  //Third: Rotation about the Z-axis of a quantity (phi)
  Rotate3DZ(phi_pos+180*TMath::DegToRad());
  //
}

void UNISFaziaQuartet::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(x_angle*TMath::RadToDeg());
}

void UNISFaziaQuartet::Rotate3DY(Double_t y_angle)
{
  fDetectorMatrix->RotateY(y_angle*TMath::RadToDeg());
}

void UNISFaziaQuartet::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISFaziaQuartet::Translate3D(Double_t x, Double_t y, Double_t z)
{  
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(x,y,z));
}

// destructor
UNISFaziaQuartet::~UNISFaziaQuartet()
{
  for(Int_t i=0; i<TRowColumn; i++)  
  {
    delete[] TCenters[i];
    delete[] TCentersXprime[i];
    delete[] TCentersYprime[i];
  }  
  delete [] TCenters;
  delete [] TCentersXprime;
  delete [] TCentersYprime;
}

// returns 1 if the particle is inside the telescope, 0 if not.
// this function sets also the impact point XY coordinates
Int_t UNISFaziaQuartet::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
  if(TMath::Abs(TImpactX)<=TTelescopeTrue_semi && TMath::Abs(TImpactY)<=TTelescopeTrue_semi)
  {
    return 1;
  }
  return 0;
}

// 
// Returns an absolute number identifying the pixel fired according to the following scheme:
// back*num_pads + front
// If the particle is not inside the active area -> return value = -1
Int_t UNISFaziaQuartet::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  Int_t i,j;
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  /*If not the particle is inside the telescope*/

  // Impact point coordinates with respect to the Top Left corner
  // These quantities are always positive within the surface of the detector
  TImpactXprime=TTopLeftXCorner-TImpactX;
  TImpactYprime=TImpactY+TTopLeftYCorner;
  // Matrix i,j indexes inside the telescope
  // i = pad front (vertical)
  // j = pad back (horizontal)
  i=Int_t(TImpactYprime/TPadTrue_width);
  j=Int_t(TImpactXprime/TPadTrue_width);
  
  /*check if the particle is inside the effective area*/
  if(fabs(TImpactXprime-TCentersXprime[i][j])<TPadEffective_semi && fabs(TImpactYprime-TCentersYprime[i][j])<TPadEffective_semi) 
  {
    return j*TRowColumn+i;
  }
  /*particle not inside the effective area*/
  return -1;
}

//draws the telescope on the X-Y plane
void UNISFaziaQuartet::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  return; 
}

//3D drawing function
void UNISFaziaQuartet::Draw3D(Option_t * draw_opt) const
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

// returns the pointer to a TGraph object that contains all the pads centers
TGraph* UNISFaziaQuartet::GetGraphObject()
{
  Double_t x[TPads_number],y[TPads_number];
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      x[i*TRowColumn+j]=TCenters[i][j].Y(); 
      y[i*TRowColumn+j]=TCenters[i][j].X();
    }
  }
  TGraph * GraphCenters = new TGraph(TPads_number,x,y);
  GraphCenters->SetMarkerStyle(20);
  GraphCenters->SetMarkerColor(kBlue);
  GraphCenters->SetMarkerSize(0.07);
  return GraphCenters;
}

// returns the TVector3 of the telescope center
TVector3 UNISFaziaQuartet::GetDetectorCenter()
{
  return TTelescopeCenter; 
}

// returns the TVector3 of the pixel center identified by a given pad front and back in the lab reference frame
TVector3 UNISFaziaQuartet::GetPadCenter(int pad)
{
  const int row=pad/TRowColumn;
  const int column=pad%TRowColumn;
  return TCenters[row][column]; 
}

TVector3 UNISFaziaQuartet::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TImpactX*TXversor+TImpactY*TYversor+TTelescopeCenter;
}

#ifdef GRAPHICAL_DEBUG
void UNISFaziaQuartet::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  Draw();

  if(IsInside(theta_inc,phi_inc,x0,y0,z0))
  {
    TGraph * impact_point = new TGraph();
    TVector3 TheImpactPoint = TImpactY*TYversor+TImpactX*TXversor+TTelescopeCenter;
    impact_point->SetPoint(0,TheImpactPoint.Y(),TheImpactPoint.X());
    impact_point->SetMarkerColor(kRed);
    impact_point->Draw("* SAME");
  }
}
#endif
