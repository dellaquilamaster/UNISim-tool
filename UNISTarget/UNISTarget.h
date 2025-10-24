#ifndef UNISTARGET_H
#define UNISTARGET_H

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


//Possible future
//set TARGET_NAME "stack_of_MgO_w_backing" //e.g. [MgO_1000ug][Al_25ug][MgO_150ug][Al_25ug]
//DEFINED in TARGET_STACK_H
//set TARGET_TILT_X 45.0 //degrees around x axis, z-axis is beam
//set TARGET_TILT_Y 15.0 //degrees around y axis, z-axis is beam

//DEFINED in TARGET_H; from beam towards detectors //both ug and um MUST be given
//add material "MgO" -thickness_ug=1000 -thickness_um=... -active=yes/no
//add material "Al" -thickness_ug=25 -thickness_um=... -active=yes/no
//add material "MgO" -thickness_ug=150 -thickness_um=... -active=yes/no
//add material "Al" -thickness_ug=25 -thickness_um=... -active=yes/no
//

typedef struct Target_
{
  std::string material;
  bool IsActive; //active/passive
  double thickness_um;
  double thickness_ug;
} Target;

class UNISTarget {
public :
  UNISTarget();
  UNISTarget(const char* name, const char* material, bool isActive,
             double thickness_um, double thickness_ug);  
  ~UNISTarget();
  
  void SetName(const char *);
  void SetMaterial(const char*);
  void SetThickness_um(double thickness_um);
  void SetThickness_ug(double thickness_ug);
  void SetIsActive(bool isActive);

  void Print() const;
  const Target& GetTarget() const;
  const std::string GetName() const;
  std::string GetTargetMaterial() const;
  double GetTargetThickness_um() const;
  double GetTargetThickness_ug() const;
  bool TargetIsActive() const;

  bool operator==(const UNISTarget& other) const {
    return GetName() == other.GetName() &&
           GetTargetMaterial() == other.GetTargetMaterial() &&
           GetTargetThickness_ug() == other.GetTargetThickness_ug() &&
           GetTargetThickness_um() == other.GetTargetThickness_um() &&
           TargetIsActive() == other.TargetIsActive();
  }

private :
  std::string fName;
  Target target;
};

#endif
