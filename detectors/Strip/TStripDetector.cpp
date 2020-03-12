#include <TVector3.h>
#include <TMath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "TStripDetector.h"

// standard constructor
TStripDetector::TStripDetector(Double_t distance, Double_t theta_pos, Double_t phi_pos, Int_t N_Strips, 
				   Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TTiltAngle(-9999)
{
  TPixels_number=N_Strips*N_Strips;
  TStrips_number=N_Strips;
  TPixelEffective_width=strip_width-inter_width;
  TInter_width=inter_width; 
  TFrame_width=frame_width; 
  TDeadLayer=dead_layer; 
  TPixelTrue_width=strip_width;
  TPixelEffective_semi=0.5*TPixelEffective_width;
  TPixelTrue_semi=0.5*TPixelTrue_width;
  TTelescopeEffective_semi=TPixelTrue_semi*TStrips_number;
  TTelescopeTrue_semi=TTelescopeEffective_semi+TDeadLayer+TFrame_width;
  TStrip_hit=(Int_t*)new Int_t[2];
  TTrueImpactPoint.SetXYZ(0.,0.,0.);
  TTelescopeImpactPoint.SetXYZ(0.,0.,0.);
  /*Telescope's corners*/
  TTopLeftCorner.SetXYZ(TTelescopeTrue_semi,-TTelescopeTrue_semi,distance);
  TTopRightCorner.SetXYZ(TTelescopeTrue_semi,TTelescopeTrue_semi,distance);
  TBottomLeftCorner.SetXYZ(-TTelescopeTrue_semi,-TTelescopeTrue_semi,distance);
  TTopLeftXCorner= TTelescopeEffective_semi;
  TTopLeftYCorner= TTelescopeEffective_semi;
  /* vectors allocations*/
  TCenters       =new TVector3*[TStrips_number];
  TCentersXprime =new Double_t*[TStrips_number];
  TCentersYprime =new Double_t*[TStrips_number];
  for(Int_t i = 0; i<TStrips_number; i++)
  {
    TCenters[i]       = new TVector3[TStrips_number]; 
    TCentersXprime[i] = new Double_t[TStrips_number]; 
    TCentersYprime[i] = new Double_t[TStrips_number]; 
  }

  /*Setting the position of the telescope center*/
  TTelescopeCenter.SetXYZ(0.,0.,distance);
  /*Rotation of the telescope center to the input position*/
  TTelescopeCenter.RotateX(theta_pos);
  TTelescopeCenter.RotateZ(phi_pos+180*TMath::DegToRad());
  Ta=TTelescopeCenter.X()/TTelescopeCenter.Mag();
  Tb=TTelescopeCenter.Y()/TTelescopeCenter.Mag();
  Tc=TTelescopeCenter.Z()/TTelescopeCenter.Mag();
  Td=-distance;
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
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
  /*************************************************************************************/
  // placement of the entire cluster of pixels on a plane orthogonal to z-axis at the input distance
  // here i represents the strip front and j represents the strip back
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].SetXYZ((TStrips_number-(2*j+1))*TPixelTrue_semi,-(TStrips_number-(2*i+1))*TPixelTrue_semi,distance);
    }
  }
  // first rotation around the z axis
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].RotateZ(-phi_pos-180*TMath::DegToRad());
    }
  }
  // rotation of the entire custer of pixels to the inpunt (theta,phi) position
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].RotateX(theta_pos);
      TCenters[i][j].RotateZ(phi_pos+180*TMath::DegToRad());
    }
  }
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
  /**************************************************************************************/
}

// constructor with arbitrary tilt angle, center position (X0, Y0, Z0)
// tilt_X = tilt angle with respect to the X-axis (vertical)
TStripDetector::TStripDetector(Double_t X0, Double_t Y0, Double_t Z0, Double_t tilt_X, 
                   Int_t N_Strips, Double_t strip_width, Double_t inter_width, Double_t frame_width, Double_t dead_layer, Option_t *opt) :
TXlabversor(1,0,0),
TYlabversor(0,1,0),
TZlabversor(0,0,1),
TTiltAngle(tilt_X),
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
  //Placing the entire cluster centered at (0,0,0) and perpendicular to the beam axis
  /*Telescope's corners*/
  TTopLeftCorner.SetXYZ(TTelescopeTrue_semi,-TTelescopeTrue_semi,0.);
  TTopRightCorner.SetXYZ(TTelescopeTrue_semi,TTelescopeTrue_semi,0.);
  TBottomLeftCorner.SetXYZ(-TTelescopeTrue_semi,-TTelescopeTrue_semi,0.);
  TTopLeftXCorner= TTelescopeEffective_semi;
  TTopLeftYCorner= TTelescopeEffective_semi;
  /* vectors allocations*/
  TCenters       =new TVector3*[TStrips_number];
  TCentersXprime =new Double_t*[TStrips_number];
  TCentersYprime =new Double_t*[TStrips_number];
  for(Int_t i = 0; i<TStrips_number; i++)
  {
    TCenters[i]       = new TVector3[TStrips_number]; 
    TCentersXprime[i] = new Double_t[TStrips_number]; 
    TCentersYprime[i] = new Double_t[TStrips_number]; 
  }
  //
  // placement of the entire cluster of pixels on a plane orthogonal to z-axis and centered at (0,0,0)
  // here i represents the strip front and j represents the strip back
  // Setting the position of the telescope center first
  TTelescopeCenter.SetXYZ(0.,0.,0.);
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j].SetXYZ((TStrips_number-(2*j+1))*TPixelTrue_semi,-(TStrips_number-(2*i+1))*TPixelTrue_semi,0.);
    }
  }
  //
  //Rotation of the detector to the input tilt angle
  //
  RotateX(TTiltAngle);
  //
  //Translation of the detector to the input position
  //
  TTelescopeCenter=TTelescopeCenter+X0*TXlabversor+Y0*TYlabversor+Z0*TZlabversor; //Detector center
  TTopLeftCorner+=X0*TXlabversor+Y0*TYlabversor+Z0*TZlabversor;
  TTopRightCorner+=X0*TXlabversor+Y0*TYlabversor+Z0*TZlabversor;
  TBottomLeftCorner+=X0*TXlabversor+Y0*TYlabversor+Z0*TZlabversor;
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {
      TCenters[i][j]=TCenters[i][j]+X0*TXlabversor+Y0*TYlabversor+Z0*TZlabversor; //i,j strip center
    }
  }
  //
  
  //
  //Calculation of the detector plane equation
  TVector3 OriginalDirection(0,0,1);
  OriginalDirection.RotateX(TTiltAngle);
  Ta=OriginalDirection.X()/OriginalDirection.Mag();
  Tb=OriginalDirection.Y()/OriginalDirection.Mag();
  Tc=OriginalDirection.Z()/OriginalDirection.Mag();
  Td=-(Ta*TTelescopeCenter.X()+Tb*TTelescopeCenter.Y()+Tc*TTelescopeCenter.Z());
  //
}

void TStripDetector::RotateX(Double_t x_angle)
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
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
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
  
  return;
}

void TStripDetector::RotateY(Double_t y_angle)
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
  /*calculation of the X Y versors*/
  TXversor=(TTopLeftCorner-TBottomLeftCorner);
  TXversor*=(1./TXversor.Mag());
  TYversor=(TTopRightCorner-TTopLeftCorner);
  TYversor*=(1./TYversor.Mag());    
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
  
  return;
}

// destructor
TStripDetector::~TStripDetector()
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
Int_t TStripDetector::IsInside(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Int_t TStripDetector::GetPixel(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Int_t TStripDetector::GetStripFront(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Int_t TStripDetector::GetStripBack(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
Double_t TStripDetector::GetThetaPixel(Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -100;
  return GetPixelCenter(numfront,numback).Theta();
}

//returns -100 if the particle is not inside the telescope
Double_t TStripDetector::GetPhiPixel(Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -100;
  return GetPixelCenter(numfront,numback).Phi();
}

//returns -1 if the particle is not inside the telescope in the other cases returns 1 and write theta and phi detected 
//in the input memory addresses
Int_t TStripDetector::GetThetaPhiPixel(Double_t * ptheta_det, Double_t * pphi_det, Int_t numfront, Int_t numback)
{
  if  (numfront<0 || numback<0) return -1;
  *ptheta_det=TStripDetector::GetThetaPixel(numfront, numback);
  *pphi_det=TStripDetector::GetPhiPixel(numfront, numback);
  return 1;
}

//draws the telescope on the X-Y plane
void TStripDetector::Draw(Option_t * draw_opt, double Xmin, double Xmax, double Ymin, double Ymax) const
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
void TStripDetector::Draw3D(Option_t * draw_opt) const
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
  
  //Generating detector frame
  TEveGeoShape * DetectorFrame = new TEveGeoShape("DetectorFrame");
  DetectorFrame->SetShape(new TGeoPara(TTelescopeTrue_semi, TTelescopeTrue_semi, 0.05, 0, 0, 0));
  DetectorFrame->SetMainColor(kYellow+2);
  //
  //Rotation (translation) of the frame to the input position
  if(TTiltAngle==-9999) {
    //First constructor used, the detector is facing the target perpendicularly
    double theta_angle=TTelescopeCenter.Theta();
    double phi_angle=TTelescopeCenter.Phi()-TMath::Pi()/2.;
    DetectorFrame->RefMainTrans().Move3PF(0., 0., TTelescopeCenter.Mag());
    DetectorFrame->RefMainTrans().RotatePF(1, 2, -phi_angle-180*TMath::DegToRad());
    DetectorFrame->RefMainTrans().RotatePF(3, 1, theta_angle);
    DetectorFrame->RefMainTrans().RotatePF(1, 2, phi_angle+180*TMath::DegToRad());
  } else {
    //Second constructor used, the detector has a tilt angle with respect to the X (vertical) axis and an absolute translaction
    DetectorFrame->RefMainTrans().RotatePF(3, 1, TTiltAngle);
    DetectorFrame->RefMainTrans().Move3PF(-TTelescopeCenter.Y(), TTelescopeCenter.X(), TTelescopeCenter.Z());
  }
  gEve->AddElement(DetectorFrame);
  //
  
  //
  //Generating pixels
  TEveGeoShape * Pixel[TStrips_number][TStrips_number];
  for(Int_t i=0; i<TStrips_number; i++)
  {
    for(Int_t j=0; j<TStrips_number; j++)
    {      
      Pixel[i][j]  = new TEveGeoShape(Form("Pixel_%02d_%02d",i,j));
      Pixel[i][j]->SetShape(new TGeoPara(TPixelEffective_semi, TPixelEffective_semi, 0.1, 0, 0, 0));
      Pixel[i][j]->SetMainColor(kGray);
      Pixel[i][j]->RefMainTrans().Move3PF((TStrips_number-(2*i+1))*TPixelTrue_semi,(TStrips_number-(2*j+1))*TPixelTrue_semi,0);
      //
      //Rotation (translation) to the final position
      if(TTiltAngle==-9999) {
        //First constructor used, the detector is facing the target perpendicularly
        double theta_angle=TTelescopeCenter.Theta();
        double phi_angle=TTelescopeCenter.Phi()-TMath::Pi()/2.;
        Pixel[i][j]->RefMainTrans().Move3PF(0., 0., TTelescopeCenter.Mag());
        Pixel[i][j]->RefMainTrans().RotatePF(1, 2, -phi_angle-180*TMath::DegToRad());
        Pixel[i][j]->RefMainTrans().RotatePF(3, 1, theta_angle);
        Pixel[i][j]->RefMainTrans().RotatePF(1, 2, phi_angle+180*TMath::DegToRad());
      } else {
        //Second constructor used, the detector has a tilt angle with respect to the X (vertical) axis and an absolute translaction
        Pixel[i][j]->RefMainTrans().RotatePF(3, 1, TTiltAngle);
        Pixel[i][j]->RefMainTrans().Move3PF(-TTelescopeCenter.Y(), TTelescopeCenter.X(), TTelescopeCenter.Z());
      }
      //
      DetectorFrame->AddElement(Pixel[i][j]);
    }
  }
  //
  
  // Drawing3D  
  gEve->Redraw3D(kTRUE);
  //
}

// returns the pointer to a TGraph object that contains all the pads centers
TGraph* TStripDetector::GetGraphObject()
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
TVector3 TStripDetector::GetDetectorCenter()
{
  return TTelescopeCenter; 
}

// returns the TVector3 of the pixel center identified by a given strip front and back in the lab reference frame
TVector3 TStripDetector::GetPixelCenter(int stripf, int stripb)
{
  return TCenters[stripf][stripb]; 
}

TVector3 TStripDetector::GetImpactPointLab(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
{
  if(!IsInside(theta_inc, phi_inc, x0, y0, z0)) return TVector3(0,0,0);
    
  return TImpactX*TXversor+TImpactY*TYversor+TTelescopeCenter;
}

#ifdef GRAPHICAL_DEBUG
void TStripDetector::ShowImpactPoint(Double_t theta_inc, Double_t phi_inc, Double_t x0, Double_t y0, Double_t z0)
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
