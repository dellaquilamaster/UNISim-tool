#include "RelativisticKinematics.h"

#include <math.h>

//_______________________________________________________________
RelativisticKinematics::RelativisticKinematics()
{}

//_______________________________________________________________
RelativisticKinematics::~RelativisticKinematics()
{}

//_______________________________________________________________
double RelativisticKinematics::GetEnergyTwoBodyEjectile(double polar_angle, double E_Projectile, double M_Projectile, double M_Target, double EStar_Projectile, double EStar_Target)
{
  double UMA=931.494; /*MeV/c^2*/
  double M_Projectile_MeV=M_Projectile*UMA;
  double M_Target_MeV=M_Target*UMA;
  double M_Ejectile_MeV=M_Projectile_MeV+EStar_Projectile;
  double M_Recoil_MeV=M_Target_MeV+EStar_Target;
  //Energia cinetica del fascio incidente
  double TKE_Projectile=E_Projectile*M_Projectile;
  //Impulso dei proiettili incidenti (modulo)
  double P_Inc=sqrt(pow(TKE_Projectile,2)+2*M_Projectile_MeV*TKE_Projectile);
  //Energia totale della particella incidente
  double E_Inc_tot=sqrt(P_Inc*P_Inc+pow(M_Projectile_MeV,2));
  //Invariant mass
  double S = pow(M_Projectile_MeV+M_Target_MeV,2)+2*M_Target_MeV*TKE_Projectile;
  //Momentum of each particle in the 0 momentum frame PRIOR
  double Pcm_prior=sqrt((pow(S-pow(M_Target_MeV,2)-pow(M_Projectile_MeV,2),2)-4*pow(M_Target_MeV,2)*pow(M_Projectile_MeV,2))/(4*S));
  //Rapity of the 0 momentum frame PRIOR
  double chi = log((Pcm_prior+sqrt(pow(M_Target_MeV,2)+pow(Pcm_prior,2)))/M_Target_MeV);
  //Momentum of each particle in the 0 momentum frame POST
  double Pcm_post=sqrt((pow(S-pow(M_Ejectile_MeV,2)-pow(M_Recoil_MeV,2),2)-4*pow(M_Ejectile_MeV,2)*pow(M_Recoil_MeV,2))/(4*S));
  //Momentum of the ejectile in the lab frame
  double P_ejectile=(sqrt(pow(M_Ejectile_MeV,2)+pow(Pcm_post,2))*cos(polar_angle)*sinh(chi)+cosh(chi)*sqrt(pow(Pcm_post,2)-pow(M_Ejectile_MeV,2)*pow(sin(polar_angle),2)*pow(sinh(chi),2)))/(1+pow(sin(polar_angle),2)*pow(sinh(chi),2));
  return sqrt(pow(P_ejectile,2)+pow(M_Ejectile_MeV,2))-M_Ejectile_MeV;
}


//_______________________________________________________________
double RelativisticKinematics::GetEnergyTwoBodyRecoil(double polar_angle, double E_Projectile, double M_Projectile, double M_Target, double EStar_Projectile, double EStar_Target)
{  
  double UMA=931.494; /*MeV/c^2*/
  double M_Projectile_MeV=M_Projectile*UMA;
  double M_Target_MeV=M_Target*UMA;
  double M_Ejectile_MeV=M_Projectile_MeV+EStar_Projectile;
  double M_Recoil_MeV=M_Target_MeV+EStar_Target;
  //Energia cinetica del fascio incidente
  double TKE_Projectile=E_Projectile*M_Projectile;
  //Impulso dei proiettili incidenti (modulo)
  double P_Inc=sqrt(pow(TKE_Projectile,2)+2*M_Projectile_MeV*TKE_Projectile);
  //Energia totale della particella incidente
  double E_Inc_tot=sqrt(P_Inc*P_Inc+pow(M_Projectile_MeV,2));
  //Invariant mass
  double S = pow(M_Projectile_MeV+M_Target_MeV,2)+2*M_Target_MeV*TKE_Projectile;
  //Momentum of each particle in the 0 momentum frame PRIOR
  double Pcm_prior=sqrt((pow(S-pow(M_Target_MeV,2)-pow(M_Projectile_MeV,2),2)-4*pow(M_Target_MeV,2)*pow(M_Projectile_MeV,2))/(4*S));
  //Rapity of the 0 momentum frame PRIOR
  double chi = log((Pcm_prior+sqrt(pow(M_Target_MeV,2)+pow(Pcm_prior,2)))/M_Target_MeV);
  //Momentum of each particle in the 0 momentum frame POST
  double Pcm_post=sqrt((pow(S-pow(M_Ejectile_MeV,2)-pow(M_Recoil_MeV,2),2)-4*pow(M_Ejectile_MeV,2)*pow(M_Recoil_MeV,2))/(4*S));
  //Momentum of the recoil in the lab frame
  double P_recoil=(sqrt(pow(M_Recoil_MeV,2)+pow(Pcm_post,2))*cos(polar_angle)*sinh(chi)+cosh(chi)*sqrt(pow(Pcm_post,2)-pow(M_Recoil_MeV,2)*pow(sin(polar_angle),2)*pow(sinh(chi),2)))/(1+pow(sin(polar_angle),2)*pow(sinh(chi),2));
  
  return sqrt(pow(P_recoil,2)+pow(M_Recoil_MeV,2))-M_Recoil_MeV;
}

//_______________________________________________________________
double RelativisticKinematics::GetEnergyTwoBodyEjectileReaction(double polar_angle, double E_Projectile, double M_Projectile, double M_Target, double M_Projectile_prime, double M_Target_prime, double EStar_Projectile, double EStar_Target)
{
  double UMA=931.494; /*MeV/c^2*/
  double M_Projectile_MeV=M_Projectile*UMA;
  double M_Target_MeV=M_Target*UMA;
  double M_Ejectile_MeV=M_Projectile_prime*UMA+EStar_Projectile;
  double M_Recoil_MeV=M_Target_prime*UMA+EStar_Target;
  //Energia cinetica del fascio incidente
  double TKE_Projectile=E_Projectile*M_Projectile;
  //Impulso dei proiettili incidenti (modulo)
  double P_Inc=sqrt(pow(TKE_Projectile,2)+2*M_Projectile_MeV*TKE_Projectile);
  //Energia totale della particella incidente
  double E_Inc_tot=sqrt(P_Inc*P_Inc+pow(M_Projectile_MeV,2));
  //Invariant mass
  double S = pow(M_Projectile_MeV+M_Target_MeV,2)+2*M_Target_MeV*TKE_Projectile;
  //Momentum of each particle in the 0 momentum frame PRIOR
  double Pcm_prior=sqrt((pow(S-pow(M_Target_MeV,2)-pow(M_Projectile_MeV,2),2)-4*pow(M_Target_MeV,2)*pow(M_Projectile_MeV,2))/(4*S));
  //Rapity of the 0 momentum frame PRIOR
  double chi = log((Pcm_prior+sqrt(pow(M_Target_MeV,2)+pow(Pcm_prior,2)))/M_Target_MeV);
  //Momentum of each particle in the 0 momentum frame POST
  double Pcm_post=sqrt((pow(S-pow(M_Ejectile_MeV,2)-pow(M_Recoil_MeV,2),2)-4*pow(M_Ejectile_MeV,2)*pow(M_Recoil_MeV,2))/(4*S));
  //Momentum of the ejectile in the lab frame
  double P_ejectile=(sqrt(pow(M_Ejectile_MeV,2)+pow(Pcm_post,2))*cos(polar_angle)*sinh(chi)+cosh(chi)*sqrt(pow(Pcm_post,2)-pow(M_Ejectile_MeV,2)*pow(sin(polar_angle),2)*pow(sinh(chi),2)))/(1+pow(sin(polar_angle),2)*pow(sinh(chi),2));
  return sqrt(pow(P_ejectile,2)+pow(M_Ejectile_MeV,2))-M_Ejectile_MeV;
}

//_______________________________________________________________
double RelativisticKinematics::GetEnergyTwoBodyRecoilReaction(double polar_angle, double E_Projectile, double M_Projectile, double M_Target, double M_Projectile_prime, double M_Target_prime, double EStar_Projectile, double EStar_Target)
{  
  double UMA=931.494; /*MeV/c^2*/
  double M_Projectile_MeV=M_Projectile*UMA;
  double M_Target_MeV=M_Target*UMA;
  double M_Ejectile_MeV=M_Projectile_prime*UMA+EStar_Projectile;
  double M_Recoil_MeV=M_Target_prime*UMA+EStar_Target;
  //Energia cinetica del fascio incidente
  double TKE_Projectile=E_Projectile*M_Projectile;
  //Impulso dei proiettili incidenti (modulo)
  double P_Inc=sqrt(pow(TKE_Projectile,2)+2*M_Projectile_MeV*TKE_Projectile);
  //Energia totale della particella incidente
  double E_Inc_tot=sqrt(P_Inc*P_Inc+pow(M_Projectile_MeV,2));
  //Invariant mass
  double S = pow(M_Projectile_MeV+M_Target_MeV,2)+2*M_Target_MeV*TKE_Projectile;
  //Momentum of each particle in the 0 momentum frame PRIOR
  double Pcm_prior=sqrt((pow(S-pow(M_Target_MeV,2)-pow(M_Projectile_MeV,2),2)-4*pow(M_Target_MeV,2)*pow(M_Projectile_MeV,2))/(4*S));
  //Rapity of the 0 momentum frame PRIOR
  double chi = log((Pcm_prior+sqrt(pow(M_Target_MeV,2)+pow(Pcm_prior,2)))/M_Target_MeV);
  //Momentum of each particle in the 0 momentum frame POST
  double Pcm_post=sqrt((pow(S-pow(M_Ejectile_MeV,2)-pow(M_Recoil_MeV,2),2)-4*pow(M_Ejectile_MeV,2)*pow(M_Recoil_MeV,2))/(4*S));
  //Momentum of the recoil in the lab frame
  double P_recoil=(sqrt(pow(M_Recoil_MeV,2)+pow(Pcm_post,2))*cos(polar_angle)*sinh(chi)+cosh(chi)*sqrt(pow(Pcm_post,2)-pow(M_Recoil_MeV,2)*pow(sin(polar_angle),2)*pow(sinh(chi),2)))/(1+pow(sin(polar_angle),2)*pow(sinh(chi),2));
  
  return sqrt(pow(P_recoil,2)+pow(M_Recoil_MeV,2))-M_Recoil_MeV;
}


