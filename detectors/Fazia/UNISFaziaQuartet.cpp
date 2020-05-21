#include "UNISFaziaQuartet.h"

// standard constructor
UNISFaziaQuartet::UNISFaziaQuartet(Double_t displacement, Double_t theta_pos, Double_t phi_pos, Double_t pad_width, Double_t frame_width) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TPads_number(4),
TRowColumn(2),
TNominalDistance(100),
TNominalTheta(theta_pos),
TNominalPhi(phi_pos),
TDisplacement(displacement),
TPadEffective_width(pad_width),
TPadEffective_semi(0.5*TPadEffective_width),
TFrame_width(frame_width),
TPadTrue_width(TPadEffective_width+TFrame_width),
TPadTrue_semi(0.5*TPadTrue_width),
TTelescopeTrue_semi(TPadTrue_semi*TRowColumn),
TTrueImpactPoint(0.,0.,0.),
TTelescopeImpactPoint(0.,0.,0.)
{
  /*Telescope's corners*/
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
  /*Setting the position of the telescope center*/
  TTelescopeCenter.SetXYZ(0.,0.,TNominalDistance);
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateX(theta_pos);
  TTelescopeCenter.RotateZ(phi_pos+180*TMath::DegToRad());
  //Moving the center to the displaced position
  TTelescopeCenter+=TVector3(0,0,TDisplacement);
  //
  Ta=TTelescopeCenter.X()/TTelescopeCenter.Mag();
  Tb=TTelescopeCenter.Y()/TTelescopeCenter.Mag();
  Tc=TTelescopeCenter.Z()/TTelescopeCenter.Mag();
  Td=-TNominalDistance;
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateZ(-phi_pos-180*TMath::DegToRad()); 
  TTopRightCorner  .RotateZ(-phi_pos-180*TMath::DegToRad()); 
  TBottomLeftCorner.RotateZ(-phi_pos-180*TMath::DegToRad());
  TTopLeftCorner   .RotateX(theta_pos);
  TTopRightCorner  .RotateX(theta_pos);
  TBottomLeftCorner.RotateX(theta_pos);
  TTopLeftCorner   .RotateZ(phi_pos+180*TMath::DegToRad());
  TTopRightCorner  .RotateZ(phi_pos+180*TMath::DegToRad());
  TBottomLeftCorner.RotateZ(phi_pos+180*TMath::DegToRad());
  //Translation to the corners to the displaced position
  TTopLeftCorner   +=TVector3(0,0,TDisplacement);
  TTopRightCorner  +=TVector3(0,0,TDisplacement);
  TBottomLeftCorner+=TVector3(0,0,TDisplacement);
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
  /*************************************************************************************/
  // placement of the entire cluster of pixels on a plane orthogonal to z-axis at the input TNominalDistance
  // here i represents the pad front and j represents the pad back
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].SetXYZ((TRowColumn-(2*j+1))*TPadTrue_semi,-(TRowColumn-(2*i+1))*TPadTrue_semi,TNominalDistance);
    }
  }
  // first rotation around the z axis
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateZ(-phi_pos-180*TMath::DegToRad());
    }
  }
  // rotation of the entire custer of pixels to the inpunt (theta,phi) position
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateX(theta_pos);
      TCenters[i][j].RotateZ(phi_pos+180*TMath::DegToRad());
    }
  }
  // Translation of the entire cluster of pads to the displaced position
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j]+=TVector3(0,0,TDisplacement);
    }
  }
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
  /**************************************************************************************/
}

void UNISFaziaQuartet::RotateX(Double_t x_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateX(x_angle);
  //
  
  /*Rotation of the pad's centers*/
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateX(x_angle);
    }
  }
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateX(x_angle); 
  TTopRightCorner  .RotateX(x_angle); 
  TBottomLeftCorner.RotateX(x_angle);
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
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
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateX(x_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  return;
}

void UNISFaziaQuartet::RotateY(Double_t y_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateY(y_angle);
  //
  /*Rotation of the pad's centers*/
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {
      TCenters[i][j].RotateY(y_angle);
    }
  }
  /*rotation of the corners to the input position*/
  TTopLeftCorner   .RotateY(y_angle); 
  TTopRightCorner  .RotateY(y_angle); 
  TBottomLeftCorner.RotateY(y_angle);
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
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
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(Ta,Tb,Tc);
  OriginalDirection.RotateY(y_angle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
  
  return;
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
  TEveGeoShape * DetectorQuartet = new TEveGeoShape("DetectorQuartet");
  gEve->AddElement(DetectorQuartet);
  //
  
  //
  //Generating pads
  TEveGeoShape * Pad[TRowColumn][TRowColumn];
  TEveGeoShape * TopFrame[TRowColumn][TRowColumn];
  TEveGeoShape * BottomFrame[TRowColumn][TRowColumn];
  TEveGeoShape * LeftFrame[TRowColumn][TRowColumn];
  TEveGeoShape * RightFrame[TRowColumn][TRowColumn];
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      Pad[i][j]  = new TEveGeoShape(Form("Pad_%02d_%02d",i,j));
      Pad[i][j]->SetShape(new TGeoPara(TPadEffective_semi, TPadEffective_semi, 0.1, 0, 0, 0));
      Pad[i][j]->SetMainColor(kGray);
      Pad[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      TopFrame[i][j]  = new TEveGeoShape(Form("TopFrame_%02d_%02d",i,j));
      TopFrame[i][j]->SetShape(new TGeoPara(TPadTrue_semi,TFrame_width/2., 0.2, 0, 0, 0));
      TopFrame[i][j]->SetMainColor(kYellow+2);
      TopFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,0);
      BottomFrame[i][j]  = new TEveGeoShape(Form("BottomFrame_%02d_%02d",i,j));
      BottomFrame[i][j]->SetShape(new TGeoPara(TPadTrue_semi,TFrame_width/2., 0.2, 0, 0, 0));
      BottomFrame[i][j]->SetMainColor(kYellow+2);
      BottomFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,0);
      LeftFrame[i][j]  = new TEveGeoShape(Form("LeftFrame_%02d_%02d",i,j));
      LeftFrame[i][j]->SetShape(new TGeoPara(TFrame_width/2.,TPadTrue_semi, 0.2, 0, 0, 0));
      LeftFrame[i][j]->SetMainColor(kYellow+2);
      LeftFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi-TPadEffective_semi-TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      RightFrame[i][j]  = new TEveGeoShape(Form("RightFrame_%02d_%02d",i,j));
      RightFrame[i][j]->SetShape(new TGeoPara(TFrame_width/2.,TPadTrue_semi, 0.2, 0, 0, 0));
      RightFrame[i][j]->SetMainColor(kYellow+2);
      RightFrame[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi+TPadEffective_semi+TFrame_width/2.,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      //
      //Rotation (translation) to the final position
      Pad[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      Pad[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      Pad[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      Pad[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      Pad[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      TopFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      TopFrame[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      TopFrame[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      TopFrame[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      TopFrame[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      BottomFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      BottomFrame[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      BottomFrame[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      BottomFrame[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      BottomFrame[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      LeftFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      LeftFrame[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      LeftFrame[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      LeftFrame[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      LeftFrame[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      RightFrame[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance);
      RightFrame[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      RightFrame[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      RightFrame[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      RightFrame[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      //
      DetectorQuartet->AddElement(Pad[i][j]);
      DetectorQuartet->AddElement(TopFrame[i][j]);
      DetectorQuartet->AddElement(BottomFrame[i][j]);
      DetectorQuartet->AddElement(LeftFrame[i][j]);
      DetectorQuartet->AddElement(RightFrame[i][j]);
    }
  }
  //Drawing decoratively a block of CsIs
  TEveGeoShape * CsICrystal[TRowColumn][TRowColumn];
  const double CsICrystalLength=10.; // cm
  for(Int_t i=0; i<TRowColumn; i++)
  {
    for(Int_t j=0; j<TRowColumn; j++)
    {      
      CsICrystal[i][j]  = new TEveGeoShape(Form("Pad_%02d_%02d",i,j));
      CsICrystal[i][j]->SetShape(new TGeoPara(TPadEffective_semi, TPadEffective_semi, CsICrystalLength, 0, 0, 0));
      CsICrystal[i][j]->SetMainColor(kBlue-4);
      CsICrystal[i][j]->RefMainTrans().Move3PF((TRowColumn-(2*i+1))*TPadTrue_semi,(TRowColumn-(2*j+1))*TPadTrue_semi,0);
      //
      //Rotation (translation) to the final position
      CsICrystal[i][j]->RefMainTrans().Move3PF(0., 0., TNominalDistance+CsICrystalLength+0.1);
      CsICrystal[i][j]->RefMainTrans().RotatePF(1, 2, -TNominalPhi-180*TMath::DegToRad());
      CsICrystal[i][j]->RefMainTrans().RotatePF(3, 1, TNominalTheta);
      CsICrystal[i][j]->RefMainTrans().RotatePF(1, 2, TNominalPhi+180*TMath::DegToRad());
      CsICrystal[i][j]->RefMainTrans().Move3PF(0., 0., TDisplacement);
      //
      DetectorQuartet->AddElement(CsICrystal[i][j]);
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
TVector3 UNISFaziaQuartet::GetPadCenter(int padf, int padb)
{
  return TCenters[padf][padb]; 
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
