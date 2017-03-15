// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIContainedROIdIContainedROIDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "ContainedROIConfig.h"
#include "ContainedROI.h"
#include "FlashROIMatching.h"

// Header files passed via #pragma extra_include

namespace larlitecv {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *larlitecv_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("larlitecv", 0 /*version*/, "ContainedROIConfig.h", 7,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &larlitecv_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *larlitecv_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace {
  void TriggerDictionaryInitialization_libLArliteCVContainedROI_Impl() {
    static const char* headers[] = {
"ContainedROIConfig.h",
"ContainedROI.h",
"FlashROIMatching.h",
0
    };
    static const char* includePaths[] = {
"/home/barnchri/CV_Directories_Test/LArCV/build/include",
"/usr/local/include",
"/usr/include/python2.7",
"/usr/include/x86_64-linux-gnu/python2.7",
"/usr/local/lib/python2.7/dist-packages/numpy/core/include",
"/home/barnchri/CV_Directories_Test/LArCV/build/include",
"/usr/local/include",
"/usr/include/python2.7",
"/usr/include/x86_64-linux-gnu/python2.7",
"/usr/local/lib/python2.7/dist-packages/numpy/core/include/../app",
"/home/barnchri/CV_Directories_Test/larlite/core",
"/home/barnchri/CV_Directories_Test/larlite/core/../UserDev/BasicTool",
"/home/barnchri/CV_Directories_Test/larlite/core/../UserDev/SelectionTool",
"/home/barnchri/CV_Directories_Test/larlitecv/build/include",
"/home/barnchri/CV_Directories_Test/LArCV/app/ann_1.1.2/include",
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/app/ContainedROI/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVContainedROI dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVContainedROI dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "ContainedROIConfig.h"
#include "ContainedROI.h"
#include "FlashROIMatching.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVContainedROI",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVContainedROI_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVContainedROI_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVContainedROI() {
  TriggerDictionaryInitialization_libLArliteCVContainedROI_Impl();
}
