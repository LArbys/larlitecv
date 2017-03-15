// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIGapChsdIGapChsDict

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
#include "EmptyChannelAlgo.h"
#include "GapChProcessor.h"

// Header files passed via #pragma extra_include

namespace larlitecv {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *larlitecv_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("larlitecv", 0 /*version*/, "GapChs/EmptyChannelAlgo.h", 15,
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

namespace ROOT {
   static TClass *larlitecvcLcLEmptyChannelAlgo_Dictionary();
   static void larlitecvcLcLEmptyChannelAlgo_TClassManip(TClass*);
   static void *new_larlitecvcLcLEmptyChannelAlgo(void *p = 0);
   static void *newArray_larlitecvcLcLEmptyChannelAlgo(Long_t size, void *p);
   static void delete_larlitecvcLcLEmptyChannelAlgo(void *p);
   static void deleteArray_larlitecvcLcLEmptyChannelAlgo(void *p);
   static void destruct_larlitecvcLcLEmptyChannelAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::EmptyChannelAlgo*)
   {
      ::larlitecv::EmptyChannelAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::EmptyChannelAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::EmptyChannelAlgo", "GapChs/EmptyChannelAlgo.h", 17,
                  typeid(::larlitecv::EmptyChannelAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLEmptyChannelAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::EmptyChannelAlgo) );
      instance.SetNew(&new_larlitecvcLcLEmptyChannelAlgo);
      instance.SetNewArray(&newArray_larlitecvcLcLEmptyChannelAlgo);
      instance.SetDelete(&delete_larlitecvcLcLEmptyChannelAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLEmptyChannelAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLEmptyChannelAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::EmptyChannelAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::EmptyChannelAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::EmptyChannelAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLEmptyChannelAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::EmptyChannelAlgo*)0x0)->GetClass();
      larlitecvcLcLEmptyChannelAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLEmptyChannelAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLGapChProcessor_Dictionary();
   static void larlitecvcLcLGapChProcessor_TClassManip(TClass*);
   static void *new_larlitecvcLcLGapChProcessor(void *p = 0);
   static void *newArray_larlitecvcLcLGapChProcessor(Long_t size, void *p);
   static void delete_larlitecvcLcLGapChProcessor(void *p);
   static void deleteArray_larlitecvcLcLGapChProcessor(void *p);
   static void destruct_larlitecvcLcLGapChProcessor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::GapChProcessor*)
   {
      ::larlitecv::GapChProcessor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::GapChProcessor));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::GapChProcessor", "GapChs/GapChProcessor.h", 11,
                  typeid(::larlitecv::GapChProcessor), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLGapChProcessor_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::GapChProcessor) );
      instance.SetNew(&new_larlitecvcLcLGapChProcessor);
      instance.SetNewArray(&newArray_larlitecvcLcLGapChProcessor);
      instance.SetDelete(&delete_larlitecvcLcLGapChProcessor);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLGapChProcessor);
      instance.SetDestructor(&destruct_larlitecvcLcLGapChProcessor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::GapChProcessor*)
   {
      return GenerateInitInstanceLocal((::larlitecv::GapChProcessor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::GapChProcessor*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLGapChProcessor_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::GapChProcessor*)0x0)->GetClass();
      larlitecvcLcLGapChProcessor_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLGapChProcessor_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLEmptyChannelAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::EmptyChannelAlgo : new ::larlitecv::EmptyChannelAlgo;
   }
   static void *newArray_larlitecvcLcLEmptyChannelAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::EmptyChannelAlgo[nElements] : new ::larlitecv::EmptyChannelAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLEmptyChannelAlgo(void *p) {
      delete ((::larlitecv::EmptyChannelAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLEmptyChannelAlgo(void *p) {
      delete [] ((::larlitecv::EmptyChannelAlgo*)p);
   }
   static void destruct_larlitecvcLcLEmptyChannelAlgo(void *p) {
      typedef ::larlitecv::EmptyChannelAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::EmptyChannelAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLGapChProcessor(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::GapChProcessor : new ::larlitecv::GapChProcessor;
   }
   static void *newArray_larlitecvcLcLGapChProcessor(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::GapChProcessor[nElements] : new ::larlitecv::GapChProcessor[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLGapChProcessor(void *p) {
      delete ((::larlitecv::GapChProcessor*)p);
   }
   static void deleteArray_larlitecvcLcLGapChProcessor(void *p) {
      delete [] ((::larlitecv::GapChProcessor*)p);
   }
   static void destruct_larlitecvcLcLGapChProcessor(void *p) {
      typedef ::larlitecv::GapChProcessor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::GapChProcessor

namespace ROOT {
   static TClass *vectorlEvectorlElarcvcLcLImage2DgRsPgR_Dictionary();
   static void vectorlEvectorlElarcvcLcLImage2DgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p);
   static void destruct_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<larcv::Image2D> >*)
   {
      vector<vector<larcv::Image2D> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<larcv::Image2D> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<larcv::Image2D> >", -2, "vector", 210,
                  typeid(vector<vector<larcv::Image2D> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlElarcvcLcLImage2DgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<larcv::Image2D> >) );
      instance.SetNew(&new_vectorlEvectorlElarcvcLcLImage2DgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlElarcvcLcLImage2DgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlElarcvcLcLImage2DgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<larcv::Image2D> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<larcv::Image2D> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlElarcvcLcLImage2DgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<larcv::Image2D> >*)0x0)->GetClass();
      vectorlEvectorlElarcvcLcLImage2DgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlElarcvcLcLImage2DgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<larcv::Image2D> > : new vector<vector<larcv::Image2D> >;
   }
   static void *newArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<larcv::Image2D> >[nElements] : new vector<vector<larcv::Image2D> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p) {
      delete ((vector<vector<larcv::Image2D> >*)p);
   }
   static void deleteArray_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p) {
      delete [] ((vector<vector<larcv::Image2D> >*)p);
   }
   static void destruct_vectorlEvectorlElarcvcLcLImage2DgRsPgR(void *p) {
      typedef vector<vector<larcv::Image2D> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<larcv::Image2D> >

namespace ROOT {
   static TClass *vectorlElarcvcLcLImage2DgR_Dictionary();
   static void vectorlElarcvcLcLImage2DgR_TClassManip(TClass*);
   static void *new_vectorlElarcvcLcLImage2DgR(void *p = 0);
   static void *newArray_vectorlElarcvcLcLImage2DgR(Long_t size, void *p);
   static void delete_vectorlElarcvcLcLImage2DgR(void *p);
   static void deleteArray_vectorlElarcvcLcLImage2DgR(void *p);
   static void destruct_vectorlElarcvcLcLImage2DgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larcv::Image2D>*)
   {
      vector<larcv::Image2D> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larcv::Image2D>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larcv::Image2D>", -2, "vector", 210,
                  typeid(vector<larcv::Image2D>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarcvcLcLImage2DgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<larcv::Image2D>) );
      instance.SetNew(&new_vectorlElarcvcLcLImage2DgR);
      instance.SetNewArray(&newArray_vectorlElarcvcLcLImage2DgR);
      instance.SetDelete(&delete_vectorlElarcvcLcLImage2DgR);
      instance.SetDeleteArray(&deleteArray_vectorlElarcvcLcLImage2DgR);
      instance.SetDestructor(&destruct_vectorlElarcvcLcLImage2DgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larcv::Image2D> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larcv::Image2D>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarcvcLcLImage2DgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larcv::Image2D>*)0x0)->GetClass();
      vectorlElarcvcLcLImage2DgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarcvcLcLImage2DgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarcvcLcLImage2DgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larcv::Image2D> : new vector<larcv::Image2D>;
   }
   static void *newArray_vectorlElarcvcLcLImage2DgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larcv::Image2D>[nElements] : new vector<larcv::Image2D>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarcvcLcLImage2DgR(void *p) {
      delete ((vector<larcv::Image2D>*)p);
   }
   static void deleteArray_vectorlElarcvcLcLImage2DgR(void *p) {
      delete [] ((vector<larcv::Image2D>*)p);
   }
   static void destruct_vectorlElarcvcLcLImage2DgR(void *p) {
      typedef vector<larcv::Image2D> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larcv::Image2D>

namespace {
  void TriggerDictionaryInitialization_libLArliteCVGapChs_Impl() {
    static const char* headers[] = {
"EmptyChannelAlgo.h",
"GapChProcessor.h",
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
"/home/barnchri/CV_Directories_Test/larlitecv/build/include",
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/app/GapChs/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVGapChs dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlitecv{class __attribute__((annotate("$clingAutoload$EmptyChannelAlgo.h")))  EmptyChannelAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$GapChProcessor.h")))  GapChProcessor;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVGapChs dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "EmptyChannelAlgo.h"
#include "GapChProcessor.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlitecv::EmptyChannelAlgo", payloadCode, "@",
"larlitecv::GapChProcessor", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVGapChs",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVGapChs_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVGapChs_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVGapChs() {
  TriggerDictionaryInitialization_libLArliteCVGapChs_Impl();
}
