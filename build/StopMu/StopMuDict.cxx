// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIStopMudIStopMuDict

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
#include "SMClusterTypes.h"
#include "StopMuAlgo.h"
#include "StopMuAlgoTypes.h"
#include "StopMuClusterConfig.h"
#include "StopMuCluster.h"
#include "StopMuFilterSpacePoints.h"
#include "StopMuSkeleton.h"
#include "StopMuStart.h"
#include "StopMuTrackerConfig.h"
#include "StopMuTracker.h"

// Header files passed via #pragma extra_include

namespace larlitecv {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *larlitecv_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("larlitecv", 0 /*version*/, "SMClusterTypes.h", 12,
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
   static TClass *larlitecvcLcLStopMuStart_Dictionary();
   static void larlitecvcLcLStopMuStart_TClassManip(TClass*);
   static void *new_larlitecvcLcLStopMuStart(void *p = 0);
   static void *newArray_larlitecvcLcLStopMuStart(Long_t size, void *p);
   static void delete_larlitecvcLcLStopMuStart(void *p);
   static void deleteArray_larlitecvcLcLStopMuStart(void *p);
   static void destruct_larlitecvcLcLStopMuStart(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::StopMuStart*)
   {
      ::larlitecv::StopMuStart *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::StopMuStart));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::StopMuStart", "StopMuStart.h", 17,
                  typeid(::larlitecv::StopMuStart), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLStopMuStart_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::StopMuStart) );
      instance.SetNew(&new_larlitecvcLcLStopMuStart);
      instance.SetNewArray(&newArray_larlitecvcLcLStopMuStart);
      instance.SetDelete(&delete_larlitecvcLcLStopMuStart);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLStopMuStart);
      instance.SetDestructor(&destruct_larlitecvcLcLStopMuStart);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::StopMuStart*)
   {
      return GenerateInitInstanceLocal((::larlitecv::StopMuStart*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::StopMuStart*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLStopMuStart_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::StopMuStart*)0x0)->GetClass();
      larlitecvcLcLStopMuStart_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLStopMuStart_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLStopMuSkeleton_Dictionary();
   static void larlitecvcLcLStopMuSkeleton_TClassManip(TClass*);
   static void *new_larlitecvcLcLStopMuSkeleton(void *p = 0);
   static void *newArray_larlitecvcLcLStopMuSkeleton(Long_t size, void *p);
   static void delete_larlitecvcLcLStopMuSkeleton(void *p);
   static void deleteArray_larlitecvcLcLStopMuSkeleton(void *p);
   static void destruct_larlitecvcLcLStopMuSkeleton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::StopMuSkeleton*)
   {
      ::larlitecv::StopMuSkeleton *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::StopMuSkeleton));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::StopMuSkeleton", "StopMuSkeleton.h", 8,
                  typeid(::larlitecv::StopMuSkeleton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLStopMuSkeleton_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::StopMuSkeleton) );
      instance.SetNew(&new_larlitecvcLcLStopMuSkeleton);
      instance.SetNewArray(&newArray_larlitecvcLcLStopMuSkeleton);
      instance.SetDelete(&delete_larlitecvcLcLStopMuSkeleton);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLStopMuSkeleton);
      instance.SetDestructor(&destruct_larlitecvcLcLStopMuSkeleton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::StopMuSkeleton*)
   {
      return GenerateInitInstanceLocal((::larlitecv::StopMuSkeleton*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::StopMuSkeleton*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLStopMuSkeleton_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::StopMuSkeleton*)0x0)->GetClass();
      larlitecvcLcLStopMuSkeleton_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLStopMuSkeleton_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLStopMuAlgo_Dictionary();
   static void larlitecvcLcLStopMuAlgo_TClassManip(TClass*);
   static void *new_larlitecvcLcLStopMuAlgo(void *p = 0);
   static void *newArray_larlitecvcLcLStopMuAlgo(Long_t size, void *p);
   static void delete_larlitecvcLcLStopMuAlgo(void *p);
   static void deleteArray_larlitecvcLcLStopMuAlgo(void *p);
   static void destruct_larlitecvcLcLStopMuAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::StopMuAlgo*)
   {
      ::larlitecv::StopMuAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::StopMuAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::StopMuAlgo", "StopMuAlgo.h", 26,
                  typeid(::larlitecv::StopMuAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLStopMuAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::StopMuAlgo) );
      instance.SetNew(&new_larlitecvcLcLStopMuAlgo);
      instance.SetNewArray(&newArray_larlitecvcLcLStopMuAlgo);
      instance.SetDelete(&delete_larlitecvcLcLStopMuAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLStopMuAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLStopMuAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::StopMuAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::StopMuAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::StopMuAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLStopMuAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::StopMuAlgo*)0x0)->GetClass();
      larlitecvcLcLStopMuAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLStopMuAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLStep3D_Dictionary();
   static void larlitecvcLcLStep3D_TClassManip(TClass*);
   static void *new_larlitecvcLcLStep3D(void *p = 0);
   static void *newArray_larlitecvcLcLStep3D(Long_t size, void *p);
   static void delete_larlitecvcLcLStep3D(void *p);
   static void deleteArray_larlitecvcLcLStep3D(void *p);
   static void destruct_larlitecvcLcLStep3D(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::Step3D*)
   {
      ::larlitecv::Step3D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::Step3D));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::Step3D", "StopMuTracker.h", 32,
                  typeid(::larlitecv::Step3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLStep3D_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::Step3D) );
      instance.SetNew(&new_larlitecvcLcLStep3D);
      instance.SetNewArray(&newArray_larlitecvcLcLStep3D);
      instance.SetDelete(&delete_larlitecvcLcLStep3D);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLStep3D);
      instance.SetDestructor(&destruct_larlitecvcLcLStep3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::Step3D*)
   {
      return GenerateInitInstanceLocal((::larlitecv::Step3D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::Step3D*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLStep3D_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::Step3D*)0x0)->GetClass();
      larlitecvcLcLStep3D_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLStep3D_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLStopMuTracker_Dictionary();
   static void larlitecvcLcLStopMuTracker_TClassManip(TClass*);
   static void delete_larlitecvcLcLStopMuTracker(void *p);
   static void deleteArray_larlitecvcLcLStopMuTracker(void *p);
   static void destruct_larlitecvcLcLStopMuTracker(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::StopMuTracker*)
   {
      ::larlitecv::StopMuTracker *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::StopMuTracker));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::StopMuTracker", "StopMuTracker.h", 166,
                  typeid(::larlitecv::StopMuTracker), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLStopMuTracker_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::StopMuTracker) );
      instance.SetDelete(&delete_larlitecvcLcLStopMuTracker);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLStopMuTracker);
      instance.SetDestructor(&destruct_larlitecvcLcLStopMuTracker);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::StopMuTracker*)
   {
      return GenerateInitInstanceLocal((::larlitecv::StopMuTracker*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::StopMuTracker*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLStopMuTracker_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::StopMuTracker*)0x0)->GetClass();
      larlitecvcLcLStopMuTracker_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLStopMuTracker_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLStopMuStart(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuStart : new ::larlitecv::StopMuStart;
   }
   static void *newArray_larlitecvcLcLStopMuStart(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuStart[nElements] : new ::larlitecv::StopMuStart[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLStopMuStart(void *p) {
      delete ((::larlitecv::StopMuStart*)p);
   }
   static void deleteArray_larlitecvcLcLStopMuStart(void *p) {
      delete [] ((::larlitecv::StopMuStart*)p);
   }
   static void destruct_larlitecvcLcLStopMuStart(void *p) {
      typedef ::larlitecv::StopMuStart current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::StopMuStart

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLStopMuSkeleton(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuSkeleton : new ::larlitecv::StopMuSkeleton;
   }
   static void *newArray_larlitecvcLcLStopMuSkeleton(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuSkeleton[nElements] : new ::larlitecv::StopMuSkeleton[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLStopMuSkeleton(void *p) {
      delete ((::larlitecv::StopMuSkeleton*)p);
   }
   static void deleteArray_larlitecvcLcLStopMuSkeleton(void *p) {
      delete [] ((::larlitecv::StopMuSkeleton*)p);
   }
   static void destruct_larlitecvcLcLStopMuSkeleton(void *p) {
      typedef ::larlitecv::StopMuSkeleton current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::StopMuSkeleton

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLStopMuAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuAlgo : new ::larlitecv::StopMuAlgo;
   }
   static void *newArray_larlitecvcLcLStopMuAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::StopMuAlgo[nElements] : new ::larlitecv::StopMuAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLStopMuAlgo(void *p) {
      delete ((::larlitecv::StopMuAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLStopMuAlgo(void *p) {
      delete [] ((::larlitecv::StopMuAlgo*)p);
   }
   static void destruct_larlitecvcLcLStopMuAlgo(void *p) {
      typedef ::larlitecv::StopMuAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::StopMuAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLStep3D(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::Step3D : new ::larlitecv::Step3D;
   }
   static void *newArray_larlitecvcLcLStep3D(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::Step3D[nElements] : new ::larlitecv::Step3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLStep3D(void *p) {
      delete ((::larlitecv::Step3D*)p);
   }
   static void deleteArray_larlitecvcLcLStep3D(void *p) {
      delete [] ((::larlitecv::Step3D*)p);
   }
   static void destruct_larlitecvcLcLStep3D(void *p) {
      typedef ::larlitecv::Step3D current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::Step3D

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLStopMuTracker(void *p) {
      delete ((::larlitecv::StopMuTracker*)p);
   }
   static void deleteArray_larlitecvcLcLStopMuTracker(void *p) {
      delete [] ((::larlitecv::StopMuTracker*)p);
   }
   static void destruct_larlitecvcLcLStopMuTracker(void *p) {
      typedef ::larlitecv::StopMuTracker current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::StopMuTracker

namespace ROOT {
   static TClass *vectorlEvectorlEvectorlEdoublegRsPgRsPgR_Dictionary();
   static void vectorlEvectorlEvectorlEdoublegRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p);
   static void destruct_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<vector<double> > >*)
   {
      vector<vector<vector<double> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<vector<double> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<vector<double> > >", -2, "vector", 210,
                  typeid(vector<vector<vector<double> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEvectorlEdoublegRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<vector<double> > >) );
      instance.SetNew(&new_vectorlEvectorlEvectorlEdoublegRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEvectorlEdoublegRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEvectorlEdoublegRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<vector<double> > > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<vector<double> > >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEvectorlEdoublegRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<vector<double> > >*)0x0)->GetClass();
      vectorlEvectorlEvectorlEdoublegRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEvectorlEdoublegRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<double> > > : new vector<vector<vector<double> > >;
   }
   static void *newArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<double> > >[nElements] : new vector<vector<vector<double> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p) {
      delete ((vector<vector<vector<double> > >*)p);
   }
   static void deleteArray_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p) {
      delete [] ((vector<vector<vector<double> > >*)p);
   }
   static void destruct_vectorlEvectorlEvectorlEdoublegRsPgRsPgR(void *p) {
      typedef vector<vector<vector<double> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<vector<double> > >

namespace ROOT {
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary();
   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEintgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<int> >*)
   {
      vector<vector<int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<int> >", -2, "vector", 210,
                  typeid(vector<vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEintgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<int> >) );
      instance.SetNew(&new_vectorlEvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<int> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<int> >*)0x0)->GetClass();
      vectorlEvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> > : new vector<vector<int> >;
   }
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> >[nElements] : new vector<vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEintgRsPgR(void *p) {
      delete ((vector<vector<int> >*)p);
   }
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p) {
      delete [] ((vector<vector<int> >*)p);
   }
   static void destruct_vectorlEvectorlEintgRsPgR(void *p) {
      typedef vector<vector<int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<int> >

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary();
   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double> >*)
   {
      vector<vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double> >", -2, "vector", 210,
                  typeid(vector<vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<double> >) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<double> >*)0x0)->GetClass();
      vectorlEvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> > : new vector<vector<double> >;
   }
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> >[nElements] : new vector<vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete ((vector<vector<double> >*)p);
   }
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete [] ((vector<vector<double> >*)p);
   }
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p) {
      typedef vector<vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double> >

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

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 210,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEdbscancLcLdbscanOutputgR_Dictionary();
   static void vectorlEdbscancLcLdbscanOutputgR_TClassManip(TClass*);
   static void *new_vectorlEdbscancLcLdbscanOutputgR(void *p = 0);
   static void *newArray_vectorlEdbscancLcLdbscanOutputgR(Long_t size, void *p);
   static void delete_vectorlEdbscancLcLdbscanOutputgR(void *p);
   static void deleteArray_vectorlEdbscancLcLdbscanOutputgR(void *p);
   static void destruct_vectorlEdbscancLcLdbscanOutputgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<dbscan::dbscanOutput>*)
   {
      vector<dbscan::dbscanOutput> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<dbscan::dbscanOutput>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<dbscan::dbscanOutput>", -2, "vector", 210,
                  typeid(vector<dbscan::dbscanOutput>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdbscancLcLdbscanOutputgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<dbscan::dbscanOutput>) );
      instance.SetNew(&new_vectorlEdbscancLcLdbscanOutputgR);
      instance.SetNewArray(&newArray_vectorlEdbscancLcLdbscanOutputgR);
      instance.SetDelete(&delete_vectorlEdbscancLcLdbscanOutputgR);
      instance.SetDeleteArray(&deleteArray_vectorlEdbscancLcLdbscanOutputgR);
      instance.SetDestructor(&destruct_vectorlEdbscancLcLdbscanOutputgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<dbscan::dbscanOutput> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<dbscan::dbscanOutput>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdbscancLcLdbscanOutputgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<dbscan::dbscanOutput>*)0x0)->GetClass();
      vectorlEdbscancLcLdbscanOutputgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdbscancLcLdbscanOutputgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdbscancLcLdbscanOutputgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dbscan::dbscanOutput> : new vector<dbscan::dbscanOutput>;
   }
   static void *newArray_vectorlEdbscancLcLdbscanOutputgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dbscan::dbscanOutput>[nElements] : new vector<dbscan::dbscanOutput>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdbscancLcLdbscanOutputgR(void *p) {
      delete ((vector<dbscan::dbscanOutput>*)p);
   }
   static void deleteArray_vectorlEdbscancLcLdbscanOutputgR(void *p) {
      delete [] ((vector<dbscan::dbscanOutput>*)p);
   }
   static void destruct_vectorlEdbscancLcLdbscanOutputgR(void *p) {
      typedef vector<dbscan::dbscanOutput> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<dbscan::dbscanOutput>

namespace {
  void TriggerDictionaryInitialization_libLArliteCVStopMu_Impl() {
    static const char* headers[] = {
"SMClusterTypes.h",
"StopMuAlgo.h",
"StopMuAlgoTypes.h",
"StopMuClusterConfig.h",
"StopMuCluster.h",
"StopMuFilterSpacePoints.h",
"StopMuSkeleton.h",
"StopMuStart.h",
"StopMuTrackerConfig.h",
"StopMuTracker.h",
0
    };
    static const char* includePaths[] = {
"./",
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
"/home/barnchri/CV_Directories_Test/larlite/core/../UserDev/BasicTool/GeoAlgo",
"/home/barnchri/CV_Directories_Test/larlitecv/build/include",
"/home/barnchri/CV_Directories_Test/LArCV/app/ann_1.1.2/include",
"/usr/local/include",
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/app/StopMu/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVStopMu dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlitecv{class __attribute__((annotate("$clingAutoload$StopMuStart.h")))  StopMuStart;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$StopMuSkeleton.h")))  StopMuSkeleton;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$StopMuAlgo.h")))  StopMuAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$StopMuTracker.h")))  Step3D;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$StopMuTracker.h")))  StopMuTracker;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVStopMu dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef USE_OPENCV
  #define USE_OPENCV 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "SMClusterTypes.h"
#include "StopMuAlgo.h"
#include "StopMuAlgoTypes.h"
#include "StopMuClusterConfig.h"
#include "StopMuCluster.h"
#include "StopMuFilterSpacePoints.h"
#include "StopMuSkeleton.h"
#include "StopMuStart.h"
#include "StopMuTrackerConfig.h"
#include "StopMuTracker.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlitecv::Step3D", payloadCode, "@",
"larlitecv::StopMuAlgo", payloadCode, "@",
"larlitecv::StopMuSkeleton", payloadCode, "@",
"larlitecv::StopMuStart", payloadCode, "@",
"larlitecv::StopMuTracker", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVStopMu",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVStopMu_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVStopMu_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVStopMu() {
  TriggerDictionaryInitialization_libLArliteCVStopMu_Impl();
}
