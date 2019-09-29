// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME UNISRootEventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/UNISRootEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_UNISRootEvent(void *p = 0);
   static void *newArray_UNISRootEvent(Long_t size, void *p);
   static void delete_UNISRootEvent(void *p);
   static void deleteArray_UNISRootEvent(void *p);
   static void destruct_UNISRootEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::UNISRootEvent*)
   {
      ::UNISRootEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::UNISRootEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("UNISRootEvent", ::UNISRootEvent::Class_Version(), "UNISRootEvent.h", 8,
                  typeid(::UNISRootEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::UNISRootEvent::Dictionary, isa_proxy, 4,
                  sizeof(::UNISRootEvent) );
      instance.SetNew(&new_UNISRootEvent);
      instance.SetNewArray(&newArray_UNISRootEvent);
      instance.SetDelete(&delete_UNISRootEvent);
      instance.SetDeleteArray(&deleteArray_UNISRootEvent);
      instance.SetDestructor(&destruct_UNISRootEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::UNISRootEvent*)
   {
      return GenerateInitInstanceLocal((::UNISRootEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::UNISRootEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr UNISRootEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *UNISRootEvent::Class_Name()
{
   return "UNISRootEvent";
}

//______________________________________________________________________________
const char *UNISRootEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UNISRootEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int UNISRootEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UNISRootEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *UNISRootEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UNISRootEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *UNISRootEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UNISRootEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void UNISRootEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class UNISRootEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(UNISRootEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(UNISRootEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_UNISRootEvent(void *p) {
      return  p ? new(p) ::UNISRootEvent : new ::UNISRootEvent;
   }
   static void *newArray_UNISRootEvent(Long_t nElements, void *p) {
      return p ? new(p) ::UNISRootEvent[nElements] : new ::UNISRootEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_UNISRootEvent(void *p) {
      delete ((::UNISRootEvent*)p);
   }
   static void deleteArray_UNISRootEvent(void *p) {
      delete [] ((::UNISRootEvent*)p);
   }
   static void destruct_UNISRootEvent(void *p) {
      typedef ::UNISRootEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::UNISRootEvent

namespace {
  void TriggerDictionaryInitialization_UNISRootEventDict_Impl() {
    static const char* headers[] = {
"include/UNISRootEvent.h",
0
    };
    static const char* includePaths[] = {
"./detectors/DetectionSetup/",
"./detectors/Strip/",
"./detectors/Lamp/",
"./LISETools/",
"./include/",
"./generator/",
"/opt/ROOT/buildroot61804/include",
"/home/daniele/Dropbox/Ricerca/Ruder_Boskovic/software/UNISim-tool/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "UNISRootEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$include/UNISRootEvent.h")))  UNISRootEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "UNISRootEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/UNISRootEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"UNISRootEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("UNISRootEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_UNISRootEventDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_UNISRootEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_UNISRootEventDict() {
  TriggerDictionaryInitialization_UNISRootEventDict_Impl();
}
