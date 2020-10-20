#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "UNISStripDetector.h"

// standard constructor
UNISStripDetector::UNISStripDetector(Double_t distance, Double_t theta_pos, Double_t phi_pos, Int_t N_Strips, 
				   Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TTrueImpactPoint(0.,0.,0.),
TTelescopeImpactPoint(0.,0.,0.),
TTiltXAngle(-9999),
TTiltYAngle(-9999),
TStrips_number(N_Strips),                                                  
TPixels_number(N_Strips*N_Strips),                                         
TPixelTrue_width(strip_width),                                             
TPixelTrue_semi(0.5*TPixelTrue_width),                                     
TInter_width(inter_width),                                                 
TFrame_width(frame_width),                                                 
TDeadLayer(dead_layer),                                                    
TPixelEffective_width(strip_width-inter_width),                            
TPixelEffective_semi(0.5*TPixelEffective_width),                           
TTelescopeEffective_semi(TPixelTrue_semi*TStrips_number),                  
TTelescopeTrue_semi(TTelescopeEffective_semi+TDeadLayer+TFrame_width),                                             
TStrip_hit((Int_t*)new Int_t[2])
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
  //vectors allocations*
  TCenters       =new TVector3*[TStrips_number];
  TCentersXprime =new Double_t*[TStrips_number];
  TCentersYprime =new Double_t*[TStrips_number];
  //
  for(Int_t i = 0; i<TStrips_number; i++)
  {
    TCenters[i]       = new TVector3[TStrips_number]; 
    TCentersXprime[i] = new Double_t[TStrips_number]; 
    TCentersYprime[i] = new Double_t[TStrips_number]; 
  }
  //
  
  //
  //calculation of the X Y versors
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());  
  //

  //
  // placement of the entire cluster of pixels on a plane orthogonal to z-axis at the input distance
  // here i represents the strip front and j represents the strip back
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].SetXYZ((TStrips_number-(2*j+1))*TPixelTrue_semi,-(TStrips_number-(2*i+1))*TPixelTrue_semi,distance);
    }
  }
  //
  
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the strip, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
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
UNISStripDetector::UNISStripDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, Double_t tilt_Y,
                   Int_t N_Strips, Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
UNISStripDetector(0., 0., 0., N_Strips, strip_width, inter_width, frame_width, dead_layer, opt)
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

void UNISStripDetector::RotateX(Double_t x_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateX(x_angle);
  //
  
  /*Rotation of the pad's centers*/
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].RotateX(x_angle);
    }
  }
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
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the strip, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
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
  
  //
  Rotate3DX(x_angle);
  //  
  
  return;
}

void UNISStripDetector::RotateY(Double_t y_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateY(y_angle);
  //
  
  /*Rotation of the pad's centers*/
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].RotateY(y_angle);
    }
  }
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
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the strip, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
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
  
  //
  Rotate3DY(y_angle);
  //  
  
  return;
}

void UNISStripDetector::RotateZ(Double_t z_angle)
{
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateZ(z_angle);
  //
  
  /*Rotation of the pad's centers*/
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].RotateZ(z_angle);
    }
  }
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
  // determination of the pixel center coordinates respect to the top right corner of the telescope
  // Here the axes are oriented towards the inner side of the strip, such as that all the coordinates are positive within the detector
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCentersXprime[i][j]=TTopLeftXCorner-(TCenters[i][j]-TTelescopeCenter).Dot(TXversor);
      TCentersYprime[i][j]=TTopLeftYCorner+(TCenters[i][j]-TTelescopeCenter).Dot(TYversor);
    }
  }
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


void UNISStripDetector::Translate(Double_t x, Double_t y, Double_t z)
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
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
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
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
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

// destructor
UNISStripDetector::~UNISStripDetector()
{
  for(Int_t i=0; i<TStrips_number; i++)  
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
Int_t UNISStripDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Int_t UNISStripDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  Int_t i,j;
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  /*If not the particle is inside the telescope*/

  // Impact point coordinates with respect to the Top Left corner
  // These quantities are always positive within the surface of the detector
  TImpactXprime=TTopLeftXCorner-TImpactX;
  TImpactYprime=TImpactY+TTopLeftYCorner;
  // Matrix i,j indexes inside the telescope
  // i = strip front (vertical)
  // j = strip back (horizontal)
  i=Int_t(TImpactYprime/TPixelTrue_width);
  j=Int_t(TImpactXprime/TPixelTrue_width);
  
  /*check if the particle is inside the effective area*/
  if(fabs(TImpactXprime-TCentersXprime[i][j])<TPixelEffective_semi && fabs(TImpactYprime-TCentersYprime[i][j])<TPixelEffective_semi) 
  {
    return j*TStrips_number+i;
  }
  /*particle not inside the effective area*/
  return -1;
}

// returns the number of the hit strip front. Returns -1 if the particle is not within the effective area.
Int_t UNISStripDetector::GetStripFront(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  //Check if the particle is inside the active detector    
    
  Int_t i,j;
  // Impact point coordinates with respect to the Top Left corner
  // These quantities are always positive within the surface of the detector
  TImpactXprime=TTopLeftXCorner-TImpactX;
  TImpactYprime=TImpactY+TTopLeftYCorner;
  
  // Matrix i,j indexes inside the telescope
  // i = strip front (vertical)
  // j = strip back (horizontal)
  i=Int_t(TImpactYprime/TPixelTrue_width);
  j=Int_t(TImpactXprime/TPixelTrue_width);
  
  /*check if the particle is inside the effective area, if so, returs the front strip number*/
  if(fabs(TImpactXprime-TCentersXprime[i][j])<TPixelEffective_semi && fabs(TImpactYprime-TCentersYprime[i][j])<TPixelEffective_semi)  return i;
  return -1;
}

// returns the number of the hit strip front. Returns -1 if the particle is not within the effective area.
Int_t UNISStripDetector::GetStripBack(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{  
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return -1;  //Check if the particle is inside the active detector     
    
  Int_t i,j;
  // Impact point coordinates with respect to the Top Left corner
  // These quantities are always positive within the surface of the detector
  TImpactXprime=TTopLeftXCorner-TImpactX;
  TImpactYprime=TImpactY+TTopLeftYCorner;
  // Matrix i,j indexes inside the telescope
  // i = strip front (vertical)
  // j = strip back (horizontal)
  i=Int_t(TImpactYprime/TPixelTrue_width);
  j=Int_t(TImpactXprime/TPixelTrue_width);
  
  /*check if the particle is inside the effective area and, if so, returs the back strip number*/
  if(fabs(TImpactXprime-TCentersXprime[i][j])<TPixelEffective_semi && fabs(TImpactYprime-TCentersYprime[i][j])<TPixelEffective_semi)  return j;
  return -1;
}

//returns -100 if the particle is not inside the telescope
Double_t UNISStripDetector::GetThetaPixel(Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -100;
  return GetPixelCenter(numfront,numback).Theta();
}

//returns -100 if the particle is not inside the telescope
Double_t UNISStripDetector::GetPhiPixel(Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -100;
  return GetPixelCenter(numfront,numback).Phi();
}

//returns -1 if the particle is not inside the telescope in the other cases returns 1 and write theta and phi detected 
//in the input memory addresses
Int_t UNISStripDetector::GetThetaPhiPixel(Double_t * ptheta_det, Double_t * pphi_det, Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -1;
  *ptheta_det=UNISStripDetector::GetThetaPixel(numfront, numback);
  *pphi_det=UNISStripDetector::GetPhiPixel(numfront, numback);
  return 1;
}

//draws the telescope on the X-Y plane
void UNISStripDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
{
  TLine * LeftOuterFrameBorder=new TLine((TCenters[0][0]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[0][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X(),(TCenters[0][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[0][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X());
  TLine * RightOuterFrameBorder=new TLine((TCenters[TStrips_number-1][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[TStrips_number-1][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X(),(TCenters[TStrips_number-1][TStrips_number-1]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[TStrips_number-1][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X());
  TLine * TopOuterFrameBorder=new TLine((TCenters[0][0]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[0][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X(),(TCenters[TStrips_number-1][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[TStrips_number-1][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X());
  TLine * BottomOuterFrameBorder=new TLine((TCenters[0][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[0][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X(),(TCenters[TStrips_number-1][TStrips_number-1]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y(),(TCenters[TStrips_number-1][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X());
  
  TLatex * FirstFront = new TLatex((TCenters[0][0]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).Y(),(TCenters[0][0]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).X(),"0");
  TLatex * LastFront = new TLatex((TCenters[TStrips_number-1][0]+TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).Y(),(TCenters[TStrips_number-1][0]+TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).X(),Form("%d",TStrips_number-1));
  TLatex * FirstBack = new TLatex((TCenters[0][0]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).Y(),(TCenters[0][0]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).X(),"0");
  TLatex * LastBack = new TLatex((TCenters[0][TStrips_number-1]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).Y(),(TCenters[0][TStrips_number-1]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi).X(),Form("%d",TStrips_number-1));

  TLine * LeftPixelBorder[TStrips_number][TStrips_number];
  TLine * RightPixelBorder[TStrips_number][TStrips_number];
  TLine * BottomPixelBorder[TStrips_number][TStrips_number];
  TLine * TopPixelBorder[TStrips_number][TStrips_number];
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TVector3 TopLeftCorner     = TCenters[i][j]-TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi;
      TVector3 TopRightCorner    = TCenters[i][j]+TYversor*TPixelEffective_semi+TXversor*TPixelEffective_semi;
      TVector3 BottomLeftCorner  = TCenters[i][j]-TYversor*TPixelEffective_semi-TXversor*TPixelEffective_semi;
      TVector3 BottomRightCorner = TCenters[i][j]+TYversor*TPixelEffective_semi-TXversor*TPixelEffective_semi;
      
      LeftPixelBorder[i][j]  =new TLine(TopLeftCorner.Y(),TopLeftCorner.X(),BottomLeftCorner.Y(),BottomLeftCorner.X()); 
      RightPixelBorder[i][j] =new TLine(TopRightCorner.Y(),TopRightCorner.X(),BottomRightCorner.Y(),BottomRightCorner.X());
      BottomPixelBorder[i][j]=new TLine(BottomLeftCorner.Y(),BottomLeftCorner.X(),BottomRightCorner.Y(),BottomRightCorner.X()); 
      TopPixelBorder[i][j]   =new TLine(TopLeftCorner.Y(),TopLeftCorner.X(),TopRightCorner.Y(),TopRightCorner.X());
    }
  }  
  
  if(strstr(draw_opt,"SAME")==0 && strstr(draw_opt,"same")==0) {
    TCanvas * c1 = new TCanvas("c1","Upstream view, (0,0) is the beamline", 600,600);
    TGraph * TheCenter = new TGraph();
    TheCenter->SetPoint(0,0,0);
    TheCenter->Draw("A*");
    if(Xmin==0 && Xmax==0 && Ymin==0 && Ymax==0) {
      TheCenter->GetXaxis()->SetLimits(-fabs((TCenters[0][0]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y())*1.2,+
                                       fabs((TCenters[TStrips_number-1][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TYversor).Y())*1.2);
      TheCenter->GetYaxis()->SetRangeUser(-fabs((TCenters[0][TStrips_number-1]-(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X())*1.2,
                                          fabs((TCenters[0][0]+(TDeadLayer+TFrame_width+TPixelTrue_semi)*TXversor).X())*1.2);
    } else {
      TheCenter->GetXaxis()->SetLimits(Xmin, Xmax);
      TheCenter->GetYaxis()->SetRangeUser(Ymin, Ymax);
    }
  }
  
  LeftOuterFrameBorder->SetLineColor(kRed);
  RightOuterFrameBorder->SetLineColor(kRed);
  TopOuterFrameBorder->SetLineColor(kRed);
  BottomOuterFrameBorder->SetLineColor(kRed);
  LeftOuterFrameBorder->Draw("SAME");
  RightOuterFrameBorder->Draw("SAME");
  TopOuterFrameBorder->Draw("SAME");
  BottomOuterFrameBorder->Draw("SAME");
  LastFront->SetTextAlign(kHAlignRight);
  FirstBack->SetTextAlign(kHAlignRight+kVAlignTop);
  LastBack->SetTextAlign(kHAlignRight+kVAlignTop);
  FirstFront->Draw("SAME");
  LastFront->Draw("SAME");
  FirstBack->Draw("SAME");
  LastBack->Draw("SAME");
  
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      LeftPixelBorder[i][j]->Draw("SAME");
      RightPixelBorder[i][j]->Draw("SAME");
      BottomPixelBorder[i][j]->Draw("SAME");
      TopPixelBorder[i][j]->Draw("SAME");
    }
  }
  return; 
}

//3D drawing function
void UNISStripDetector::Draw3D(Option_t * draw_opt) const
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

void UNISStripDetector::Generate3D()
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
  fFrame       = new TGeoVolume("pad_volume",new TGeoBBox(TTelescopeTrue_semi, TTelescopeTrue_semi, 0.05));
  fPixel       = new TGeoVolume("top_frame_volume",new TGeoBBox(TPixelEffective_semi, TPixelEffective_semi, 0.1));
  fFrame->SetLineColor(kYellow+2);
  fPixel->SetLineColor(kGray);
  //
  
  //
  //Adding frame to mother volume
  fDetector->AddNode(fFrame,0,new TGeoTranslation(0., 0., 0.));
  //
  
  //
  //Pixels
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {      
      fDetector->AddNode(fPixel,i*TStrips_number+j,new TGeoTranslation(-(TStrips_number-(2*i+1))*TPixelTrue_semi,(TStrips_number-(2*j+1))*TPixelTrue_semi,0));
    }
  }
  //
}  

void UNISStripDetector::Rotate3DX(Double_t x_angle)
{
  fDetectorMatrix->RotateX(x_angle*TMath::RadToDeg());
}

void UNISStripDetector::Rotate3DY(Double_t y_angle)
{
  fDetectorMatrix->RotateY(y_angle*TMath::RadToDeg());
}

void UNISStripDetector::Rotate3DZ(Double_t z_angle)
{
  fDetectorMatrix->RotateZ(z_angle*TMath::RadToDeg());
}

void UNISStripDetector::Translate3D(Double_t x, Double_t y, Double_t z)
{  
  fDetectorMatrix->MultiplyLeft(TGeoTranslation(x,y,z));
}

// returns the pointer to a TGraph object that contains all the pads centers
TGraph* UNISStripDetector::GetGraphObject()
{
  Double_t x[TPixels_number],y[TPixels_number];
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      x[i*TStrips_number+j]=TCenters[i][j].Y(); 
      y[i*TStrips_number+j]=TCenters[i][j].X();
    }
  }
  TGraph * GraphCenters = new TGraph(TPixels_number,x,y);
  GraphCenters->SetMarkerStyle(20);
  GraphCenters->SetMarkerColor(kBlue);
  GraphCenters->SetMarkerSize(0.07);
  return GraphCenters;
}

// returns the TVector3 of the telescope center
TVector3 UNISStripDetector::GetDetectorCenter()
{
  return TTelescopeCenter; 
}

// returns the TVector3 of the pixel center identified by a given strip front and back in the lab reference frame
TVector3 UNISStripDetector::GetPixelCenter(int stripf, int stripb)
{
  return TCenters[stripf][stripb]; 
}

TVector3 UNISStripDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TImpactX*TXversor+TImpactY*TYversor+TTelescopeCenter;
}

#ifdef GRAPHICAL_DEBUG
void UNISStripDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
