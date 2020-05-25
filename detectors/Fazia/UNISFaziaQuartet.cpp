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
TPadTrue_width(TPadEffective_width+TFrame_width),
TPadTrue_semi(0.5*TPadTrue_width),
TTelescopeTrue_semi(TPadTrue_semi*TRowColumn),
TTrueImpactPoint(0.,0.,0.),
TTelescopeImpactPoint(0.,0.,0.)
{
  //
  //Placing the detector at 0 degrees and distance equale to the nominal one (100 cm)
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
  Generate3D(0., 0.);
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
  TVector3 TranslationVector (y, x, z);
  
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

void UNISFaziaQuartet::Generate3D(Double_t theta_pos, Double_t phi_pos)
{
  //Creating the Graphics Manager
  TEveManager::Create();
  //
  
  //
  DetectorQuartet = new TEveGeoShape("DetectorQuartet");
  //
  
  //
  fPad = new TEveGeoShape **[TRowColumn];
  fTopFrame = new TEveGeoShape **[TRowColumn];
  fBottomFrame = new TEveGeoShape **[TRowColumn];
  fLeftFrame = new TEveGeoShape **[TRowColumn];
  fRightFrame = new TEveGeoShape **[TRowColumn];
  //
  
  //
  for(Int_t i=0; i<TRowColumn; i++)
  {
    fPad[i] = new TEveGeoShape *[TRowColumn];
    fTopFrame[i] = new TEveGeoShape *[TRowColumn];
    fBottomFrame[i] = new TEveGeoShape *[TRowColumn];
    fLeftFrame[i] = new TEveGeoShape *[TRowColumn];
    fRightFrame[i] = new TEveGeoShape *[TRowColumn];
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      fPad[i][j]  = new TEveGeoShape(Form("fPad_%02d_%02d",i,j));
      fPad[i][j]->SetShape(new TGeoPara(TPadEffective_semi, TPadEffective_semi, 0.1, 0, 0, 0));
      fPad[i][j]->SetMainColor(kGray);
      fPad[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      fTopFrame[i][j]  = new TEveGeoShape(Form("fTopFrame_%02d_%02d",i,j));
      fTopFrame[i][j]->SetShape(new TGeoPara(TPadTrue_semi,TFrame_width/2., 0.2, 0, 0, 0));
      fTopFrame[i][j]->SetMainColor(kYellow+2);
      fTopFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,0);
      fBottomFrame[i][j]  = new TEveGeoShape(Form("fBottomFrame_%02d_%02d",i,j));
      fBottomFrame[i][j]->SetShape(new TGeoPara(TPadTrue_semi,TFrame_width/2., 0.2, 0, 0, 0));
      fBottomFrame[i][j]->SetMainColor(kYellow+2);
      fBottomFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,0);
      fLeftFrame[i][j]  = new TEveGeoShape(Form("fLeftFrame_%02d_%02d",i,j));
      fLeftFrame[i][j]->SetShape(new TGeoPara(TFrame_width/2.,TPadTrue_semi, 0.2, 0, 0, 0));
      fLeftFrame[i][j]->SetMainColor(kYellow+2);
      fLeftFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      fRightFrame[i][j]  = new TEveGeoShape(Form("fRightFrame_%02d_%02d",i,j));
      fRightFrame[i][j]->SetShape(new TGeoPara(TFrame_width/2.,TPadTrue_semi, 0.2, 0, 0, 0));
      fRightFrame[i][j]->SetMainColor(kYellow+2);
      fRightFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      //
      //Placing at 0 degrees
      fPad[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      fTopFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      fBottomFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      fLeftFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      fRightFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      //
    }
  }
  //Drawing decoratively a block of CsIs
  fCsICrystal = new TEveGeoShape **[TRowColumn];
  const double fCsICrystalLength=10.; // cm
  for(Int_t i=0; i<TRowColumn; i++)
  {
    fCsICrystal[i] = new TEveGeoShape *[TRowColumn];
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      fCsICrystal[i][j]  = new TEveGeoShape(Form("fPad_%02d_%02d",i,j));
      fCsICrystal[i][j]->SetShape(new TGeoPara(TPadEffective_semi, TPadEffective_semi, fCsICrystalLength/2., 0, 0, 0));
      fCsICrystal[i][j]->SetMainColor(kBlue-4);
      fCsICrystal[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      //
      //Placing at 0 degrees
      fCsICrystal[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance+fCsICrystalLength/2.+0.1);
      //
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
  //
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      //
      //Rotation about X axis
      fPad[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      fTopFrame[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      fBottomFrame[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      fLeftFrame[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      fRightFrame[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      //
      fCsICrystal[i][j]->RefMainTrans().RotatePF(3, 1, x_angle);
      //
    }
  }
  //
}

void UNISFaziaQuartet::Rotate3DY(Double_t y_angle)
{
  //
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      //
      //Rotation about X axis
      fPad[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      fTopFrame[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      fBottomFrame[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      fLeftFrame[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      fRightFrame[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      //
      fCsICrystal[i][j]->RefMainTrans().RotatePF(3, 2, y_angle);
      //
    }
  }
  //
}

void UNISFaziaQuartet::Rotate3DZ(Double_t z_angle)
{
  //
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      //
      //Rotation about X axis
      fPad[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      fTopFrame[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      fBottomFrame[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      fLeftFrame[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      fRightFrame[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      //
      fCsICrystal[i][j]->RefMainTrans().RotatePF(1, 2, z_angle);
      //
    }
  }
  //
}

void UNISFaziaQuartet::Translate3D(Double_t x, Double_t y, Double_t z)
{  
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      //
      //Rotation about X axis
      fPad[i][j]->RefMainTrans().Move3PF(-y, x, z);
      fTopFrame[i][j]->RefMainTrans().Move3PF(-y, x, z);
      fBottomFrame[i][j]->RefMainTrans().Move3PF(-y, x, z);
      fLeftFrame[i][j]->RefMainTrans().Move3PF(-y, x, z);
      fRightFrame[i][j]->RefMainTrans().Move3PF(-y, x, z);
      //
      fCsICrystal[i][j]->RefMainTrans().Move3PF(-y, x, z);
      //
    }
  }
  //
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
    XaxisLine->SetMainTransparency(0);
    TheAxes->AddElement(XaxisLine);
    TEveGeoShape * YaxisLine = new TEveGeoShape("YaxisLine");
    YaxisLine->SetShape(new TGeoTube(0, 0.2, 10));
    YaxisLine->SetMainColor(kRed);
    YaxisLine->RefMainTrans().RotatePF(3, 1, -90*TMath::DegToRad());
    YaxisLine->RefMainTrans().Move3PF(-10, 0, 0);
    YaxisLine->SetMainTransparency(0);
    TheAxes->AddElement(YaxisLine);
    TEveGeoShape * ZaxisLine = new TEveGeoShape("ZaxisLine");
    ZaxisLine->SetShape(new TGeoTube(0, 0.2, 10));
    ZaxisLine->SetMainColor(kBlue);
    ZaxisLine->RefMainTrans().Move3PF(0, 0, 10);
    ZaxisLine->SetMainTransparency(0);
    TheAxes->AddElement(ZaxisLine);
    TEveGeoShape * XaxisArrow = new TEveGeoShape("XaxisArrow");
    XaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    XaxisArrow->SetMainColor(kGreen+1);
    XaxisArrow->RefMainTrans().RotatePF(3, 2, 90*TMath::DegToRad());
    XaxisArrow->RefMainTrans().Move3PF(0, 20, 0);
    XaxisArrow->SetMainTransparency(0);
    TheAxes->AddElement(XaxisArrow);
    TEveGeoShape * YaxisArrow = new TEveGeoShape("YaxisArrow");
    YaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    YaxisArrow->SetMainColor(kRed);
    YaxisArrow->RefMainTrans().RotatePF(3, 1, -90*TMath::DegToRad());
    YaxisArrow->RefMainTrans().Move3PF(-20, 0, 0);
    YaxisArrow->SetMainTransparency(0);
    TheAxes->AddElement(YaxisArrow);
    TEveGeoShape * ZaxisArrow = new TEveGeoShape("ZaxisArrow");
    ZaxisArrow->SetShape(new TGeoCone(0.6, 0, 0.8, 0, 0));
    ZaxisArrow->SetMainColor(kBlue);
    ZaxisArrow->RefMainTrans().Move3PF(0, 0, 20);
    ZaxisArrow->SetMainTransparency(0);
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
  gEve->AddElement(DetectorQuartet);
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      //
      DetectorQuartet->AddElement(fPad[i][j]);
      DetectorQuartet->AddElement(fTopFrame[i][j]);
      DetectorQuartet->AddElement(fBottomFrame[i][j]);
      DetectorQuartet->AddElement(fLeftFrame[i][j]);
      DetectorQuartet->AddElement(fRightFrame[i][j]);
      DetectorQuartet->AddElement(fCsICrystal[i][j]);
      //
    }
  }
  //
  
  // Drawing3D  
  gEve->Redraw3D(kTRUE);
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
