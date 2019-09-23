#include <UNISEventGenerator.h>

//____________________________________________________
UNISEventGenerator::UNISEventGenerator() :
fRandomGen(new TRandom3 (gRandomSeed))
{}

//____________________________________________________
UNISEventGenerator::~UNISEventGenerator()
{}

//____________________________________________________
int UNISEventGenerator::LoadConfiguration(const char * file_name)
{
  return 0;
}

//____________________________________________________
std::vector<UNISIon> UNISEventGenerator::GetEvent()
{
  std::vector<UNISIon> AnEvent;
  return AnEvent;
}

//____________________________________________________
void UNISEventGenerator::SetBeam(UNISIon & TheParticle)
{
  fTheBeam=TheParticle;
}

//____________________________________________________
void UNISEventGenerator::SetTarget(UNISIon & TheParticle)
{
  fTheTarget=TheParticle;
}
