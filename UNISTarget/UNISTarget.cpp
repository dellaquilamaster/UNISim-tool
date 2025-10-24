#include "UNISTarget.h"

//_______________________________________
UNISTarget::UNISTarget() :
target{"", false, 0, 0}
{}

//_______________________________________
UNISTarget::UNISTarget(const char* name, const char* material, bool isActive,
                      double thickness_um, double thickness_ug) :
fName(name),
target{material, isActive, thickness_um, thickness_ug}
{}

//_______________________________________
UNISTarget::~UNISTarget()
{}

//_______________________________________
void UNISTarget::SetName(const char* targetName)
{
  fName = targetName;
}
  
//_______________________________________
void UNISTarget::SetMaterial(const char* material)
{
  target.material = material;
}

//_______________________________________
void UNISTarget::SetThickness_um(double thickness_um)
{
  target.thickness_um = thickness_um;
}

//_______________________________________
void UNISTarget::SetThickness_ug(double thickness_ug)
{
  target.thickness_ug = thickness_ug;
}

//_______________________________________
void UNISTarget::SetIsActive(bool isActive)
{
  target.IsActive = isActive;
}

//_______________________________________
void UNISTarget::Print() const {
  std::cout << " Target Name: " << fName
            << "\n Material: " << target.material
            << "\n Activity: " << target.IsActive
            << "\n Thickness: " << target.thickness_um << " um, "
            << target.thickness_ug << " ug/cm2" << std::endl;
}

//_______________________________________
const std::string UNISTarget::GetName() const
{
  return fName;
}

//_______________________________________
const Target& UNISTarget::GetTarget() const
{
  //printf("Target Name: %s\n", fName.c_str());
  return target;
}

//_______________________________________
std::string UNISTarget::GetTargetMaterial() const
{
  return target.material;
  //printf("Material: %s\n", target.material);
}

//_______________________________________
bool UNISTarget::TargetIsActive() const
{
  return target.IsActive;
  //printf("Specification: %s\n", target.specification);
}

//_______________________________________
double UNISTarget::GetTargetThickness_ug() const
{
  return target.thickness_ug;
  //printf("Target Thickness_ug: %f\n", target.thickness_ug);
}

//_______________________________________
double UNISTarget::GetTargetThickness_um() const
{
  return target.thickness_um;
  //printf("Target Thickness_um: %f\n", target.thickness_um);
}
