// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIThruMudIThruMuDict

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
#include "AStar3DAlgo.h"
#include "AStar3DPostProcessor.h"
#include "AStarDirAlgo.h"
#include "AStarGridAlgo.h"
#include "BezierCurve.h"
#include "BMTrackCluster2D.h"
#include "BMTrackCluster3D.h"
#include "BoundaryEndPt.h"
#include "BoundaryIntersectionAlgo.h"
#include "BoundaryMatchAlgo.h"
#include "BoundaryMatchArrays.h"
#include "BoundaryMuonTaggerAlgo.h"
#include "BoundaryMuonTaggerTypes.h"
#include "BoundarySpacePoint.h"
#include "ConfigBoundaryMuonTaggerAlgo.h"
#include "DBClusterMergerAlgo.h"
#include "EndPointFilter.h"
#include "FlashMuonTaggerAlgo.h"
#include "FlashMuonTaggerConfig.h"
#include "Linear3DChargeTagger.h"
#include "Linear3DChargeTaggerTypes.h"
#include "Linear3DFitter.h"
#include "Linear3DPostProcessor_crbversion.h"
#include "Linear3DPostProcessor.h"
#include "LineRegionTest.h"
#include "TrackTests.h"

// Header files passed via #pragma extra_include

namespace larlitecv {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *larlitecv_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("larlitecv", 0 /*version*/, "AStar3DAlgo.h", 26,
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
   static TClass *larlitecvcLcLBoundaryMatchArrays_Dictionary();
   static void larlitecvcLcLBoundaryMatchArrays_TClassManip(TClass*);
   static void delete_larlitecvcLcLBoundaryMatchArrays(void *p);
   static void deleteArray_larlitecvcLcLBoundaryMatchArrays(void *p);
   static void destruct_larlitecvcLcLBoundaryMatchArrays(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BoundaryMatchArrays*)
   {
      ::larlitecv::BoundaryMatchArrays *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BoundaryMatchArrays));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BoundaryMatchArrays", "BoundaryMatchArrays.h", 33230,
                  typeid(::larlitecv::BoundaryMatchArrays), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBoundaryMatchArrays_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BoundaryMatchArrays) );
      instance.SetDelete(&delete_larlitecvcLcLBoundaryMatchArrays);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBoundaryMatchArrays);
      instance.SetDestructor(&destruct_larlitecvcLcLBoundaryMatchArrays);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BoundaryMatchArrays*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BoundaryMatchArrays*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BoundaryMatchArrays*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBoundaryMatchArrays_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BoundaryMatchArrays*)0x0)->GetClass();
      larlitecvcLcLBoundaryMatchArrays_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBoundaryMatchArrays_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLConfigBoundaryMuonTaggerAlgo_Dictionary();
   static void larlitecvcLcLConfigBoundaryMuonTaggerAlgo_TClassManip(TClass*);
   static void *new_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p = 0);
   static void *newArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(Long_t size, void *p);
   static void delete_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p);
   static void deleteArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p);
   static void destruct_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::ConfigBoundaryMuonTaggerAlgo*)
   {
      ::larlitecv::ConfigBoundaryMuonTaggerAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::ConfigBoundaryMuonTaggerAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::ConfigBoundaryMuonTaggerAlgo", "ConfigBoundaryMuonTaggerAlgo.h", 13,
                  typeid(::larlitecv::ConfigBoundaryMuonTaggerAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLConfigBoundaryMuonTaggerAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::ConfigBoundaryMuonTaggerAlgo) );
      instance.SetNew(&new_larlitecvcLcLConfigBoundaryMuonTaggerAlgo);
      instance.SetNewArray(&newArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo);
      instance.SetDelete(&delete_larlitecvcLcLConfigBoundaryMuonTaggerAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLConfigBoundaryMuonTaggerAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::ConfigBoundaryMuonTaggerAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::ConfigBoundaryMuonTaggerAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::ConfigBoundaryMuonTaggerAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLConfigBoundaryMuonTaggerAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::ConfigBoundaryMuonTaggerAlgo*)0x0)->GetClass();
      larlitecvcLcLConfigBoundaryMuonTaggerAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLConfigBoundaryMuonTaggerAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBoundaryEndPt_Dictionary();
   static void larlitecvcLcLBoundaryEndPt_TClassManip(TClass*);
   static void *new_larlitecvcLcLBoundaryEndPt(void *p = 0);
   static void *newArray_larlitecvcLcLBoundaryEndPt(Long_t size, void *p);
   static void delete_larlitecvcLcLBoundaryEndPt(void *p);
   static void deleteArray_larlitecvcLcLBoundaryEndPt(void *p);
   static void destruct_larlitecvcLcLBoundaryEndPt(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BoundaryEndPt*)
   {
      ::larlitecv::BoundaryEndPt *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BoundaryEndPt));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BoundaryEndPt", "BoundaryEndPt.h", 9,
                  typeid(::larlitecv::BoundaryEndPt), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBoundaryEndPt_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BoundaryEndPt) );
      instance.SetNew(&new_larlitecvcLcLBoundaryEndPt);
      instance.SetNewArray(&newArray_larlitecvcLcLBoundaryEndPt);
      instance.SetDelete(&delete_larlitecvcLcLBoundaryEndPt);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBoundaryEndPt);
      instance.SetDestructor(&destruct_larlitecvcLcLBoundaryEndPt);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BoundaryEndPt*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BoundaryEndPt*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BoundaryEndPt*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBoundaryEndPt_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BoundaryEndPt*)0x0)->GetClass();
      larlitecvcLcLBoundaryEndPt_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBoundaryEndPt_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBoundarySpacePoint_Dictionary();
   static void larlitecvcLcLBoundarySpacePoint_TClassManip(TClass*);
   static void *new_larlitecvcLcLBoundarySpacePoint(void *p = 0);
   static void *newArray_larlitecvcLcLBoundarySpacePoint(Long_t size, void *p);
   static void delete_larlitecvcLcLBoundarySpacePoint(void *p);
   static void deleteArray_larlitecvcLcLBoundarySpacePoint(void *p);
   static void destruct_larlitecvcLcLBoundarySpacePoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BoundarySpacePoint*)
   {
      ::larlitecv::BoundarySpacePoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BoundarySpacePoint));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BoundarySpacePoint", "BoundarySpacePoint.h", 20,
                  typeid(::larlitecv::BoundarySpacePoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBoundarySpacePoint_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BoundarySpacePoint) );
      instance.SetNew(&new_larlitecvcLcLBoundarySpacePoint);
      instance.SetNewArray(&newArray_larlitecvcLcLBoundarySpacePoint);
      instance.SetDelete(&delete_larlitecvcLcLBoundarySpacePoint);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBoundarySpacePoint);
      instance.SetDestructor(&destruct_larlitecvcLcLBoundarySpacePoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BoundarySpacePoint*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BoundarySpacePoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BoundarySpacePoint*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBoundarySpacePoint_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BoundarySpacePoint*)0x0)->GetClass();
      larlitecvcLcLBoundarySpacePoint_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBoundarySpacePoint_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBMTrackCluster2D_Dictionary();
   static void larlitecvcLcLBMTrackCluster2D_TClassManip(TClass*);
   static void *new_larlitecvcLcLBMTrackCluster2D(void *p = 0);
   static void *newArray_larlitecvcLcLBMTrackCluster2D(Long_t size, void *p);
   static void delete_larlitecvcLcLBMTrackCluster2D(void *p);
   static void deleteArray_larlitecvcLcLBMTrackCluster2D(void *p);
   static void destruct_larlitecvcLcLBMTrackCluster2D(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BMTrackCluster2D*)
   {
      ::larlitecv::BMTrackCluster2D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BMTrackCluster2D));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BMTrackCluster2D", "BMTrackCluster2D.h", 15,
                  typeid(::larlitecv::BMTrackCluster2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBMTrackCluster2D_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BMTrackCluster2D) );
      instance.SetNew(&new_larlitecvcLcLBMTrackCluster2D);
      instance.SetNewArray(&newArray_larlitecvcLcLBMTrackCluster2D);
      instance.SetDelete(&delete_larlitecvcLcLBMTrackCluster2D);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBMTrackCluster2D);
      instance.SetDestructor(&destruct_larlitecvcLcLBMTrackCluster2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BMTrackCluster2D*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BMTrackCluster2D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BMTrackCluster2D*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBMTrackCluster2D_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BMTrackCluster2D*)0x0)->GetClass();
      larlitecvcLcLBMTrackCluster2D_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBMTrackCluster2D_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBMTrackCluster3D_Dictionary();
   static void larlitecvcLcLBMTrackCluster3D_TClassManip(TClass*);
   static void *new_larlitecvcLcLBMTrackCluster3D(void *p = 0);
   static void *newArray_larlitecvcLcLBMTrackCluster3D(Long_t size, void *p);
   static void delete_larlitecvcLcLBMTrackCluster3D(void *p);
   static void deleteArray_larlitecvcLcLBMTrackCluster3D(void *p);
   static void destruct_larlitecvcLcLBMTrackCluster3D(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BMTrackCluster3D*)
   {
      ::larlitecv::BMTrackCluster3D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BMTrackCluster3D));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BMTrackCluster3D", "BMTrackCluster3D.h", 17,
                  typeid(::larlitecv::BMTrackCluster3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBMTrackCluster3D_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BMTrackCluster3D) );
      instance.SetNew(&new_larlitecvcLcLBMTrackCluster3D);
      instance.SetNewArray(&newArray_larlitecvcLcLBMTrackCluster3D);
      instance.SetDelete(&delete_larlitecvcLcLBMTrackCluster3D);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBMTrackCluster3D);
      instance.SetDestructor(&destruct_larlitecvcLcLBMTrackCluster3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BMTrackCluster3D*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BMTrackCluster3D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BMTrackCluster3D*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBMTrackCluster3D_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BMTrackCluster3D*)0x0)->GetClass();
      larlitecvcLcLBMTrackCluster3D_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBMTrackCluster3D_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBoundaryIntersectionAlgo_Dictionary();
   static void larlitecvcLcLBoundaryIntersectionAlgo_TClassManip(TClass*);
   static void *new_larlitecvcLcLBoundaryIntersectionAlgo(void *p = 0);
   static void *newArray_larlitecvcLcLBoundaryIntersectionAlgo(Long_t size, void *p);
   static void delete_larlitecvcLcLBoundaryIntersectionAlgo(void *p);
   static void deleteArray_larlitecvcLcLBoundaryIntersectionAlgo(void *p);
   static void destruct_larlitecvcLcLBoundaryIntersectionAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BoundaryIntersectionAlgo*)
   {
      ::larlitecv::BoundaryIntersectionAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BoundaryIntersectionAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BoundaryIntersectionAlgo", "BoundaryIntersectionAlgo.h", 9,
                  typeid(::larlitecv::BoundaryIntersectionAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBoundaryIntersectionAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BoundaryIntersectionAlgo) );
      instance.SetNew(&new_larlitecvcLcLBoundaryIntersectionAlgo);
      instance.SetNewArray(&newArray_larlitecvcLcLBoundaryIntersectionAlgo);
      instance.SetDelete(&delete_larlitecvcLcLBoundaryIntersectionAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBoundaryIntersectionAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLBoundaryIntersectionAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BoundaryIntersectionAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BoundaryIntersectionAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BoundaryIntersectionAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBoundaryIntersectionAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BoundaryIntersectionAlgo*)0x0)->GetClass();
      larlitecvcLcLBoundaryIntersectionAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBoundaryIntersectionAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBoundaryMuonTaggerAlgo_Dictionary();
   static void larlitecvcLcLBoundaryMuonTaggerAlgo_TClassManip(TClass*);
   static void *new_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p = 0);
   static void *newArray_larlitecvcLcLBoundaryMuonTaggerAlgo(Long_t size, void *p);
   static void delete_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p);
   static void deleteArray_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p);
   static void destruct_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BoundaryMuonTaggerAlgo*)
   {
      ::larlitecv::BoundaryMuonTaggerAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BoundaryMuonTaggerAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BoundaryMuonTaggerAlgo", "BoundaryMuonTaggerAlgo.h", 26,
                  typeid(::larlitecv::BoundaryMuonTaggerAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBoundaryMuonTaggerAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BoundaryMuonTaggerAlgo) );
      instance.SetNew(&new_larlitecvcLcLBoundaryMuonTaggerAlgo);
      instance.SetNewArray(&newArray_larlitecvcLcLBoundaryMuonTaggerAlgo);
      instance.SetDelete(&delete_larlitecvcLcLBoundaryMuonTaggerAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBoundaryMuonTaggerAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLBoundaryMuonTaggerAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BoundaryMuonTaggerAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BoundaryMuonTaggerAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BoundaryMuonTaggerAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBoundaryMuonTaggerAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BoundaryMuonTaggerAlgo*)0x0)->GetClass();
      larlitecvcLcLBoundaryMuonTaggerAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBoundaryMuonTaggerAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLFlashMuonTaggerConfig_Dictionary();
   static void larlitecvcLcLFlashMuonTaggerConfig_TClassManip(TClass*);
   static void *new_larlitecvcLcLFlashMuonTaggerConfig(void *p = 0);
   static void *newArray_larlitecvcLcLFlashMuonTaggerConfig(Long_t size, void *p);
   static void delete_larlitecvcLcLFlashMuonTaggerConfig(void *p);
   static void deleteArray_larlitecvcLcLFlashMuonTaggerConfig(void *p);
   static void destruct_larlitecvcLcLFlashMuonTaggerConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::FlashMuonTaggerConfig*)
   {
      ::larlitecv::FlashMuonTaggerConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::FlashMuonTaggerConfig));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::FlashMuonTaggerConfig", "FlashMuonTaggerConfig.h", 11,
                  typeid(::larlitecv::FlashMuonTaggerConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLFlashMuonTaggerConfig_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::FlashMuonTaggerConfig) );
      instance.SetNew(&new_larlitecvcLcLFlashMuonTaggerConfig);
      instance.SetNewArray(&newArray_larlitecvcLcLFlashMuonTaggerConfig);
      instance.SetDelete(&delete_larlitecvcLcLFlashMuonTaggerConfig);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLFlashMuonTaggerConfig);
      instance.SetDestructor(&destruct_larlitecvcLcLFlashMuonTaggerConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::FlashMuonTaggerConfig*)
   {
      return GenerateInitInstanceLocal((::larlitecv::FlashMuonTaggerConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::FlashMuonTaggerConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLFlashMuonTaggerConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::FlashMuonTaggerConfig*)0x0)->GetClass();
      larlitecvcLcLFlashMuonTaggerConfig_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLFlashMuonTaggerConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLFlashMuonTaggerAlgo_Dictionary();
   static void larlitecvcLcLFlashMuonTaggerAlgo_TClassManip(TClass*);
   static void delete_larlitecvcLcLFlashMuonTaggerAlgo(void *p);
   static void deleteArray_larlitecvcLcLFlashMuonTaggerAlgo(void *p);
   static void destruct_larlitecvcLcLFlashMuonTaggerAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::FlashMuonTaggerAlgo*)
   {
      ::larlitecv::FlashMuonTaggerAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::FlashMuonTaggerAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::FlashMuonTaggerAlgo", "FlashMuonTaggerAlgo.h", 23,
                  typeid(::larlitecv::FlashMuonTaggerAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLFlashMuonTaggerAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::FlashMuonTaggerAlgo) );
      instance.SetDelete(&delete_larlitecvcLcLFlashMuonTaggerAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLFlashMuonTaggerAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLFlashMuonTaggerAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::FlashMuonTaggerAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::FlashMuonTaggerAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::FlashMuonTaggerAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLFlashMuonTaggerAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::FlashMuonTaggerAlgo*)0x0)->GetClass();
      larlitecvcLcLFlashMuonTaggerAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLFlashMuonTaggerAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarNode_Dictionary();
   static void larlitecvcLcLAStarNode_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStarNode(void *p = 0);
   static void *newArray_larlitecvcLcLAStarNode(Long_t size, void *p);
   static void delete_larlitecvcLcLAStarNode(void *p);
   static void deleteArray_larlitecvcLcLAStarNode(void *p);
   static void destruct_larlitecvcLcLAStarNode(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarNode*)
   {
      ::larlitecv::AStarNode *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarNode));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarNode", "AStarGridAlgo.h", 21,
                  typeid(::larlitecv::AStarNode), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarNode_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarNode) );
      instance.SetNew(&new_larlitecvcLcLAStarNode);
      instance.SetNewArray(&newArray_larlitecvcLcLAStarNode);
      instance.SetDelete(&delete_larlitecvcLcLAStarNode);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarNode);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarNode);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarNode*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarNode*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarNode*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarNode_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarNode*)0x0)->GetClass();
      larlitecvcLcLAStarNode_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarNode_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_Dictionary();
   static void priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_TClassManip(TClass*);
   static void *new_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p = 0);
   static void *newArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(Long_t size, void *p);
   static void delete_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p);
   static void deleteArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p);
   static void destruct_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >*)
   {
      ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >));
      static ::ROOT::TGenericClassInfo 
         instance("priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >", "queue", 367,
                  typeid(::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >) );
      instance.SetNew(&new_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR);
      instance.SetNewArray(&newArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR);
      instance.SetDelete(&delete_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR);
      instance.SetDeleteArray(&deleteArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR);
      instance.SetDestructor(&destruct_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR);

      ::ROOT::AddClassAlternate("priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >","priority_queue<larlitecv::AStarNode>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >*)0x0)->GetClass();
      priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarSet_Dictionary();
   static void larlitecvcLcLAStarSet_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStarSet(void *p = 0);
   static void *newArray_larlitecvcLcLAStarSet(Long_t size, void *p);
   static void delete_larlitecvcLcLAStarSet(void *p);
   static void deleteArray_larlitecvcLcLAStarSet(void *p);
   static void destruct_larlitecvcLcLAStarSet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarSet*)
   {
      ::larlitecv::AStarSet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarSet));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarSet", "AStarGridAlgo.h", 77,
                  typeid(::larlitecv::AStarSet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarSet_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarSet) );
      instance.SetNew(&new_larlitecvcLcLAStarSet);
      instance.SetNewArray(&newArray_larlitecvcLcLAStarSet);
      instance.SetDelete(&delete_larlitecvcLcLAStarSet);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarSet);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarSet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarSet*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarSet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarSet*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarSet_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarSet*)0x0)->GetClass();
      larlitecvcLcLAStarSet_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarSet_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarAlgoConfig_Dictionary();
   static void larlitecvcLcLAStarAlgoConfig_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStarAlgoConfig(void *p = 0);
   static void *newArray_larlitecvcLcLAStarAlgoConfig(Long_t size, void *p);
   static void delete_larlitecvcLcLAStarAlgoConfig(void *p);
   static void deleteArray_larlitecvcLcLAStarAlgoConfig(void *p);
   static void destruct_larlitecvcLcLAStarAlgoConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarAlgoConfig*)
   {
      ::larlitecv::AStarAlgoConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarAlgoConfig));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarAlgoConfig", "AStarGridAlgo.h", 108,
                  typeid(::larlitecv::AStarAlgoConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarAlgoConfig_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarAlgoConfig) );
      instance.SetNew(&new_larlitecvcLcLAStarAlgoConfig);
      instance.SetNewArray(&newArray_larlitecvcLcLAStarAlgoConfig);
      instance.SetDelete(&delete_larlitecvcLcLAStarAlgoConfig);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarAlgoConfig);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarAlgoConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarAlgoConfig*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarAlgoConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarAlgoConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarAlgoConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarAlgoConfig*)0x0)->GetClass();
      larlitecvcLcLAStarAlgoConfig_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarAlgoConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarGridAlgo_Dictionary();
   static void larlitecvcLcLAStarGridAlgo_TClassManip(TClass*);
   static void delete_larlitecvcLcLAStarGridAlgo(void *p);
   static void deleteArray_larlitecvcLcLAStarGridAlgo(void *p);
   static void destruct_larlitecvcLcLAStarGridAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarGridAlgo*)
   {
      ::larlitecv::AStarGridAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarGridAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarGridAlgo", "AStarGridAlgo.h", 122,
                  typeid(::larlitecv::AStarGridAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarGridAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarGridAlgo) );
      instance.SetDelete(&delete_larlitecvcLcLAStarGridAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarGridAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarGridAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarGridAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarGridAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarGridAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarGridAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarGridAlgo*)0x0)->GetClass();
      larlitecvcLcLAStarGridAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarGridAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLPointInfo_Dictionary();
   static void larlitecvcLcLPointInfo_TClassManip(TClass*);
   static void *new_larlitecvcLcLPointInfo(void *p = 0);
   static void *newArray_larlitecvcLcLPointInfo(Long_t size, void *p);
   static void delete_larlitecvcLcLPointInfo(void *p);
   static void deleteArray_larlitecvcLcLPointInfo(void *p);
   static void destruct_larlitecvcLcLPointInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::PointInfo*)
   {
      ::larlitecv::PointInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::PointInfo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::PointInfo", "Linear3DChargeTaggerTypes.h", 8,
                  typeid(::larlitecv::PointInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLPointInfo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::PointInfo) );
      instance.SetNew(&new_larlitecvcLcLPointInfo);
      instance.SetNewArray(&newArray_larlitecvcLcLPointInfo);
      instance.SetDelete(&delete_larlitecvcLcLPointInfo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLPointInfo);
      instance.SetDestructor(&destruct_larlitecvcLcLPointInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::PointInfo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::PointInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::PointInfo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLPointInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::PointInfo*)0x0)->GetClass();
      larlitecvcLcLPointInfo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLPointInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLPointInfoList_Dictionary();
   static void larlitecvcLcLPointInfoList_TClassManip(TClass*);
   static void *new_larlitecvcLcLPointInfoList(void *p = 0);
   static void *newArray_larlitecvcLcLPointInfoList(Long_t size, void *p);
   static void delete_larlitecvcLcLPointInfoList(void *p);
   static void deleteArray_larlitecvcLcLPointInfoList(void *p);
   static void destruct_larlitecvcLcLPointInfoList(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::PointInfoList*)
   {
      ::larlitecv::PointInfoList *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::PointInfoList));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::PointInfoList", "Linear3DChargeTaggerTypes.h", 37,
                  typeid(::larlitecv::PointInfoList), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLPointInfoList_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::PointInfoList) );
      instance.SetNew(&new_larlitecvcLcLPointInfoList);
      instance.SetNewArray(&newArray_larlitecvcLcLPointInfoList);
      instance.SetDelete(&delete_larlitecvcLcLPointInfoList);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLPointInfoList);
      instance.SetDestructor(&destruct_larlitecvcLcLPointInfoList);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::PointInfoList*)
   {
      return GenerateInitInstanceLocal((::larlitecv::PointInfoList*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::PointInfoList*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLPointInfoList_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::PointInfoList*)0x0)->GetClass();
      larlitecvcLcLPointInfoList_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLPointInfoList_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLLinear3DFitterConfig_Dictionary();
   static void larlitecvcLcLLinear3DFitterConfig_TClassManip(TClass*);
   static void *new_larlitecvcLcLLinear3DFitterConfig(void *p = 0);
   static void *newArray_larlitecvcLcLLinear3DFitterConfig(Long_t size, void *p);
   static void delete_larlitecvcLcLLinear3DFitterConfig(void *p);
   static void deleteArray_larlitecvcLcLLinear3DFitterConfig(void *p);
   static void destruct_larlitecvcLcLLinear3DFitterConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::Linear3DFitterConfig*)
   {
      ::larlitecv::Linear3DFitterConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::Linear3DFitterConfig));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::Linear3DFitterConfig", "Linear3DFitter.h", 18,
                  typeid(::larlitecv::Linear3DFitterConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLLinear3DFitterConfig_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::Linear3DFitterConfig) );
      instance.SetNew(&new_larlitecvcLcLLinear3DFitterConfig);
      instance.SetNewArray(&newArray_larlitecvcLcLLinear3DFitterConfig);
      instance.SetDelete(&delete_larlitecvcLcLLinear3DFitterConfig);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLLinear3DFitterConfig);
      instance.SetDestructor(&destruct_larlitecvcLcLLinear3DFitterConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::Linear3DFitterConfig*)
   {
      return GenerateInitInstanceLocal((::larlitecv::Linear3DFitterConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::Linear3DFitterConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLLinear3DFitterConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::Linear3DFitterConfig*)0x0)->GetClass();
      larlitecvcLcLLinear3DFitterConfig_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLLinear3DFitterConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLLinear3DFitter_Dictionary();
   static void larlitecvcLcLLinear3DFitter_TClassManip(TClass*);
   static void delete_larlitecvcLcLLinear3DFitter(void *p);
   static void deleteArray_larlitecvcLcLLinear3DFitter(void *p);
   static void destruct_larlitecvcLcLLinear3DFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::Linear3DFitter*)
   {
      ::larlitecv::Linear3DFitter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::Linear3DFitter));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::Linear3DFitter", "Linear3DFitter.h", 33,
                  typeid(::larlitecv::Linear3DFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLLinear3DFitter_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::Linear3DFitter) );
      instance.SetDelete(&delete_larlitecvcLcLLinear3DFitter);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLLinear3DFitter);
      instance.SetDestructor(&destruct_larlitecvcLcLLinear3DFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::Linear3DFitter*)
   {
      return GenerateInitInstanceLocal((::larlitecv::Linear3DFitter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::Linear3DFitter*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLLinear3DFitter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::Linear3DFitter*)0x0)->GetClass();
      larlitecvcLcLLinear3DFitter_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLLinear3DFitter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarDirNode_Dictionary();
   static void larlitecvcLcLAStarDirNode_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStarDirNode(void *p = 0);
   static void *newArray_larlitecvcLcLAStarDirNode(Long_t size, void *p);
   static void delete_larlitecvcLcLAStarDirNode(void *p);
   static void deleteArray_larlitecvcLcLAStarDirNode(void *p);
   static void destruct_larlitecvcLcLAStarDirNode(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarDirNode*)
   {
      ::larlitecv::AStarDirNode *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarDirNode));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarDirNode", "AStarDirAlgo.h", 24,
                  typeid(::larlitecv::AStarDirNode), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarDirNode_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarDirNode) );
      instance.SetNew(&new_larlitecvcLcLAStarDirNode);
      instance.SetNewArray(&newArray_larlitecvcLcLAStarDirNode);
      instance.SetDelete(&delete_larlitecvcLcLAStarDirNode);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarDirNode);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarDirNode);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarDirNode*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarDirNode*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarDirNode*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarDirNode_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarDirNode*)0x0)->GetClass();
      larlitecvcLcLAStarDirNode_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarDirNode_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarDirAlgoConfig_Dictionary();
   static void larlitecvcLcLAStarDirAlgoConfig_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStarDirAlgoConfig(void *p = 0);
   static void *newArray_larlitecvcLcLAStarDirAlgoConfig(Long_t size, void *p);
   static void delete_larlitecvcLcLAStarDirAlgoConfig(void *p);
   static void deleteArray_larlitecvcLcLAStarDirAlgoConfig(void *p);
   static void destruct_larlitecvcLcLAStarDirAlgoConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarDirAlgoConfig*)
   {
      ::larlitecv::AStarDirAlgoConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarDirAlgoConfig));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarDirAlgoConfig", "AStarDirAlgo.h", 169,
                  typeid(::larlitecv::AStarDirAlgoConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarDirAlgoConfig_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarDirAlgoConfig) );
      instance.SetNew(&new_larlitecvcLcLAStarDirAlgoConfig);
      instance.SetNewArray(&newArray_larlitecvcLcLAStarDirAlgoConfig);
      instance.SetDelete(&delete_larlitecvcLcLAStarDirAlgoConfig);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarDirAlgoConfig);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarDirAlgoConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarDirAlgoConfig*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarDirAlgoConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarDirAlgoConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarDirAlgoConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarDirAlgoConfig*)0x0)->GetClass();
      larlitecvcLcLAStarDirAlgoConfig_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarDirAlgoConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStarDirAlgo_Dictionary();
   static void larlitecvcLcLAStarDirAlgo_TClassManip(TClass*);
   static void delete_larlitecvcLcLAStarDirAlgo(void *p);
   static void deleteArray_larlitecvcLcLAStarDirAlgo(void *p);
   static void destruct_larlitecvcLcLAStarDirAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStarDirAlgo*)
   {
      ::larlitecv::AStarDirAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStarDirAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStarDirAlgo", "AStarDirAlgo.h", 184,
                  typeid(::larlitecv::AStarDirAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStarDirAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStarDirAlgo) );
      instance.SetDelete(&delete_larlitecvcLcLAStarDirAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStarDirAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLAStarDirAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStarDirAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStarDirAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStarDirAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStarDirAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStarDirAlgo*)0x0)->GetClass();
      larlitecvcLcLAStarDirAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStarDirAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStar3DNode_Dictionary();
   static void larlitecvcLcLAStar3DNode_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStar3DNode(void *p = 0);
   static void *newArray_larlitecvcLcLAStar3DNode(Long_t size, void *p);
   static void delete_larlitecvcLcLAStar3DNode(void *p);
   static void deleteArray_larlitecvcLcLAStar3DNode(void *p);
   static void destruct_larlitecvcLcLAStar3DNode(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStar3DNode*)
   {
      ::larlitecv::AStar3DNode *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStar3DNode));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStar3DNode", "AStar3DAlgo.h", 43,
                  typeid(::larlitecv::AStar3DNode), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStar3DNode_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStar3DNode) );
      instance.SetNew(&new_larlitecvcLcLAStar3DNode);
      instance.SetNewArray(&newArray_larlitecvcLcLAStar3DNode);
      instance.SetDelete(&delete_larlitecvcLcLAStar3DNode);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStar3DNode);
      instance.SetDestructor(&destruct_larlitecvcLcLAStar3DNode);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStar3DNode*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStar3DNode*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStar3DNode*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStar3DNode_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStar3DNode*)0x0)->GetClass();
      larlitecvcLcLAStar3DNode_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStar3DNode_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStar3DAlgoConfig_Dictionary();
   static void larlitecvcLcLAStar3DAlgoConfig_TClassManip(TClass*);
   static void *new_larlitecvcLcLAStar3DAlgoConfig(void *p = 0);
   static void *newArray_larlitecvcLcLAStar3DAlgoConfig(Long_t size, void *p);
   static void delete_larlitecvcLcLAStar3DAlgoConfig(void *p);
   static void deleteArray_larlitecvcLcLAStar3DAlgoConfig(void *p);
   static void destruct_larlitecvcLcLAStar3DAlgoConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStar3DAlgoConfig*)
   {
      ::larlitecv::AStar3DAlgoConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStar3DAlgoConfig));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStar3DAlgoConfig", "AStar3DAlgo.h", 239,
                  typeid(::larlitecv::AStar3DAlgoConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStar3DAlgoConfig_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStar3DAlgoConfig) );
      instance.SetNew(&new_larlitecvcLcLAStar3DAlgoConfig);
      instance.SetNewArray(&newArray_larlitecvcLcLAStar3DAlgoConfig);
      instance.SetDelete(&delete_larlitecvcLcLAStar3DAlgoConfig);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStar3DAlgoConfig);
      instance.SetDestructor(&destruct_larlitecvcLcLAStar3DAlgoConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStar3DAlgoConfig*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStar3DAlgoConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStar3DAlgoConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStar3DAlgoConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStar3DAlgoConfig*)0x0)->GetClass();
      larlitecvcLcLAStar3DAlgoConfig_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStar3DAlgoConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLAStar3DAlgo_Dictionary();
   static void larlitecvcLcLAStar3DAlgo_TClassManip(TClass*);
   static void delete_larlitecvcLcLAStar3DAlgo(void *p);
   static void deleteArray_larlitecvcLcLAStar3DAlgo(void *p);
   static void destruct_larlitecvcLcLAStar3DAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::AStar3DAlgo*)
   {
      ::larlitecv::AStar3DAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::AStar3DAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::AStar3DAlgo", "AStar3DAlgo.h", 265,
                  typeid(::larlitecv::AStar3DAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLAStar3DAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::AStar3DAlgo) );
      instance.SetDelete(&delete_larlitecvcLcLAStar3DAlgo);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLAStar3DAlgo);
      instance.SetDestructor(&destruct_larlitecvcLcLAStar3DAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::AStar3DAlgo*)
   {
      return GenerateInitInstanceLocal((::larlitecv::AStar3DAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::AStar3DAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLAStar3DAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::AStar3DAlgo*)0x0)->GetClass();
      larlitecvcLcLAStar3DAlgo_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLAStar3DAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLBezierCurve_Dictionary();
   static void larlitecvcLcLBezierCurve_TClassManip(TClass*);
   static void delete_larlitecvcLcLBezierCurve(void *p);
   static void deleteArray_larlitecvcLcLBezierCurve(void *p);
   static void destruct_larlitecvcLcLBezierCurve(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::BezierCurve*)
   {
      ::larlitecv::BezierCurve *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::BezierCurve));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::BezierCurve", "BezierCurve.h", 9,
                  typeid(::larlitecv::BezierCurve), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLBezierCurve_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::BezierCurve) );
      instance.SetDelete(&delete_larlitecvcLcLBezierCurve);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLBezierCurve);
      instance.SetDestructor(&destruct_larlitecvcLcLBezierCurve);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::BezierCurve*)
   {
      return GenerateInitInstanceLocal((::larlitecv::BezierCurve*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::BezierCurve*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLBezierCurve_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::BezierCurve*)0x0)->GetClass();
      larlitecvcLcLBezierCurve_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLBezierCurve_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBoundaryMatchArrays(void *p) {
      delete ((::larlitecv::BoundaryMatchArrays*)p);
   }
   static void deleteArray_larlitecvcLcLBoundaryMatchArrays(void *p) {
      delete [] ((::larlitecv::BoundaryMatchArrays*)p);
   }
   static void destruct_larlitecvcLcLBoundaryMatchArrays(void *p) {
      typedef ::larlitecv::BoundaryMatchArrays current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BoundaryMatchArrays

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::ConfigBoundaryMuonTaggerAlgo : new ::larlitecv::ConfigBoundaryMuonTaggerAlgo;
   }
   static void *newArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::ConfigBoundaryMuonTaggerAlgo[nElements] : new ::larlitecv::ConfigBoundaryMuonTaggerAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p) {
      delete ((::larlitecv::ConfigBoundaryMuonTaggerAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p) {
      delete [] ((::larlitecv::ConfigBoundaryMuonTaggerAlgo*)p);
   }
   static void destruct_larlitecvcLcLConfigBoundaryMuonTaggerAlgo(void *p) {
      typedef ::larlitecv::ConfigBoundaryMuonTaggerAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::ConfigBoundaryMuonTaggerAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBoundaryEndPt(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryEndPt : new ::larlitecv::BoundaryEndPt;
   }
   static void *newArray_larlitecvcLcLBoundaryEndPt(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryEndPt[nElements] : new ::larlitecv::BoundaryEndPt[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBoundaryEndPt(void *p) {
      delete ((::larlitecv::BoundaryEndPt*)p);
   }
   static void deleteArray_larlitecvcLcLBoundaryEndPt(void *p) {
      delete [] ((::larlitecv::BoundaryEndPt*)p);
   }
   static void destruct_larlitecvcLcLBoundaryEndPt(void *p) {
      typedef ::larlitecv::BoundaryEndPt current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BoundaryEndPt

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBoundarySpacePoint(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundarySpacePoint : new ::larlitecv::BoundarySpacePoint;
   }
   static void *newArray_larlitecvcLcLBoundarySpacePoint(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundarySpacePoint[nElements] : new ::larlitecv::BoundarySpacePoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBoundarySpacePoint(void *p) {
      delete ((::larlitecv::BoundarySpacePoint*)p);
   }
   static void deleteArray_larlitecvcLcLBoundarySpacePoint(void *p) {
      delete [] ((::larlitecv::BoundarySpacePoint*)p);
   }
   static void destruct_larlitecvcLcLBoundarySpacePoint(void *p) {
      typedef ::larlitecv::BoundarySpacePoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BoundarySpacePoint

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBMTrackCluster2D(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BMTrackCluster2D : new ::larlitecv::BMTrackCluster2D;
   }
   static void *newArray_larlitecvcLcLBMTrackCluster2D(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BMTrackCluster2D[nElements] : new ::larlitecv::BMTrackCluster2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBMTrackCluster2D(void *p) {
      delete ((::larlitecv::BMTrackCluster2D*)p);
   }
   static void deleteArray_larlitecvcLcLBMTrackCluster2D(void *p) {
      delete [] ((::larlitecv::BMTrackCluster2D*)p);
   }
   static void destruct_larlitecvcLcLBMTrackCluster2D(void *p) {
      typedef ::larlitecv::BMTrackCluster2D current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BMTrackCluster2D

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBMTrackCluster3D(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BMTrackCluster3D : new ::larlitecv::BMTrackCluster3D;
   }
   static void *newArray_larlitecvcLcLBMTrackCluster3D(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BMTrackCluster3D[nElements] : new ::larlitecv::BMTrackCluster3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBMTrackCluster3D(void *p) {
      delete ((::larlitecv::BMTrackCluster3D*)p);
   }
   static void deleteArray_larlitecvcLcLBMTrackCluster3D(void *p) {
      delete [] ((::larlitecv::BMTrackCluster3D*)p);
   }
   static void destruct_larlitecvcLcLBMTrackCluster3D(void *p) {
      typedef ::larlitecv::BMTrackCluster3D current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BMTrackCluster3D

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBoundaryIntersectionAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryIntersectionAlgo : new ::larlitecv::BoundaryIntersectionAlgo;
   }
   static void *newArray_larlitecvcLcLBoundaryIntersectionAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryIntersectionAlgo[nElements] : new ::larlitecv::BoundaryIntersectionAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBoundaryIntersectionAlgo(void *p) {
      delete ((::larlitecv::BoundaryIntersectionAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLBoundaryIntersectionAlgo(void *p) {
      delete [] ((::larlitecv::BoundaryIntersectionAlgo*)p);
   }
   static void destruct_larlitecvcLcLBoundaryIntersectionAlgo(void *p) {
      typedef ::larlitecv::BoundaryIntersectionAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BoundaryIntersectionAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryMuonTaggerAlgo : new ::larlitecv::BoundaryMuonTaggerAlgo;
   }
   static void *newArray_larlitecvcLcLBoundaryMuonTaggerAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::BoundaryMuonTaggerAlgo[nElements] : new ::larlitecv::BoundaryMuonTaggerAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p) {
      delete ((::larlitecv::BoundaryMuonTaggerAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p) {
      delete [] ((::larlitecv::BoundaryMuonTaggerAlgo*)p);
   }
   static void destruct_larlitecvcLcLBoundaryMuonTaggerAlgo(void *p) {
      typedef ::larlitecv::BoundaryMuonTaggerAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BoundaryMuonTaggerAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLFlashMuonTaggerConfig(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::FlashMuonTaggerConfig : new ::larlitecv::FlashMuonTaggerConfig;
   }
   static void *newArray_larlitecvcLcLFlashMuonTaggerConfig(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::FlashMuonTaggerConfig[nElements] : new ::larlitecv::FlashMuonTaggerConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLFlashMuonTaggerConfig(void *p) {
      delete ((::larlitecv::FlashMuonTaggerConfig*)p);
   }
   static void deleteArray_larlitecvcLcLFlashMuonTaggerConfig(void *p) {
      delete [] ((::larlitecv::FlashMuonTaggerConfig*)p);
   }
   static void destruct_larlitecvcLcLFlashMuonTaggerConfig(void *p) {
      typedef ::larlitecv::FlashMuonTaggerConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::FlashMuonTaggerConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLFlashMuonTaggerAlgo(void *p) {
      delete ((::larlitecv::FlashMuonTaggerAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLFlashMuonTaggerAlgo(void *p) {
      delete [] ((::larlitecv::FlashMuonTaggerAlgo*)p);
   }
   static void destruct_larlitecvcLcLFlashMuonTaggerAlgo(void *p) {
      typedef ::larlitecv::FlashMuonTaggerAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::FlashMuonTaggerAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStarNode(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarNode : new ::larlitecv::AStarNode;
   }
   static void *newArray_larlitecvcLcLAStarNode(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarNode[nElements] : new ::larlitecv::AStarNode[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarNode(void *p) {
      delete ((::larlitecv::AStarNode*)p);
   }
   static void deleteArray_larlitecvcLcLAStarNode(void *p) {
      delete [] ((::larlitecv::AStarNode*)p);
   }
   static void destruct_larlitecvcLcLAStarNode(void *p) {
      typedef ::larlitecv::AStarNode current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarNode

namespace ROOT {
   // Wrappers around operator new
   static void *new_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> > : new ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >;
   }
   static void *newArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >[nElements] : new ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p) {
      delete ((::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >*)p);
   }
   static void deleteArray_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p) {
      delete [] ((::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >*)p);
   }
   static void destruct_priority_queuelElarlitecvcLcLAStarNodecOvectorlElarlitecvcLcLAStarNodegRcOlesslElarlitecvcLcLAStarNodegRsPgR(void *p) {
      typedef ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::priority_queue<larlitecv::AStarNode,vector<larlitecv::AStarNode>,less<larlitecv::AStarNode> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStarSet(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarSet : new ::larlitecv::AStarSet;
   }
   static void *newArray_larlitecvcLcLAStarSet(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarSet[nElements] : new ::larlitecv::AStarSet[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarSet(void *p) {
      delete ((::larlitecv::AStarSet*)p);
   }
   static void deleteArray_larlitecvcLcLAStarSet(void *p) {
      delete [] ((::larlitecv::AStarSet*)p);
   }
   static void destruct_larlitecvcLcLAStarSet(void *p) {
      typedef ::larlitecv::AStarSet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarSet

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStarAlgoConfig(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarAlgoConfig : new ::larlitecv::AStarAlgoConfig;
   }
   static void *newArray_larlitecvcLcLAStarAlgoConfig(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarAlgoConfig[nElements] : new ::larlitecv::AStarAlgoConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarAlgoConfig(void *p) {
      delete ((::larlitecv::AStarAlgoConfig*)p);
   }
   static void deleteArray_larlitecvcLcLAStarAlgoConfig(void *p) {
      delete [] ((::larlitecv::AStarAlgoConfig*)p);
   }
   static void destruct_larlitecvcLcLAStarAlgoConfig(void *p) {
      typedef ::larlitecv::AStarAlgoConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarAlgoConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarGridAlgo(void *p) {
      delete ((::larlitecv::AStarGridAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLAStarGridAlgo(void *p) {
      delete [] ((::larlitecv::AStarGridAlgo*)p);
   }
   static void destruct_larlitecvcLcLAStarGridAlgo(void *p) {
      typedef ::larlitecv::AStarGridAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarGridAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLPointInfo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::PointInfo : new ::larlitecv::PointInfo;
   }
   static void *newArray_larlitecvcLcLPointInfo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::PointInfo[nElements] : new ::larlitecv::PointInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLPointInfo(void *p) {
      delete ((::larlitecv::PointInfo*)p);
   }
   static void deleteArray_larlitecvcLcLPointInfo(void *p) {
      delete [] ((::larlitecv::PointInfo*)p);
   }
   static void destruct_larlitecvcLcLPointInfo(void *p) {
      typedef ::larlitecv::PointInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::PointInfo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLPointInfoList(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::PointInfoList : new ::larlitecv::PointInfoList;
   }
   static void *newArray_larlitecvcLcLPointInfoList(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::PointInfoList[nElements] : new ::larlitecv::PointInfoList[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLPointInfoList(void *p) {
      delete ((::larlitecv::PointInfoList*)p);
   }
   static void deleteArray_larlitecvcLcLPointInfoList(void *p) {
      delete [] ((::larlitecv::PointInfoList*)p);
   }
   static void destruct_larlitecvcLcLPointInfoList(void *p) {
      typedef ::larlitecv::PointInfoList current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::PointInfoList

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLLinear3DFitterConfig(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::Linear3DFitterConfig : new ::larlitecv::Linear3DFitterConfig;
   }
   static void *newArray_larlitecvcLcLLinear3DFitterConfig(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::Linear3DFitterConfig[nElements] : new ::larlitecv::Linear3DFitterConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLLinear3DFitterConfig(void *p) {
      delete ((::larlitecv::Linear3DFitterConfig*)p);
   }
   static void deleteArray_larlitecvcLcLLinear3DFitterConfig(void *p) {
      delete [] ((::larlitecv::Linear3DFitterConfig*)p);
   }
   static void destruct_larlitecvcLcLLinear3DFitterConfig(void *p) {
      typedef ::larlitecv::Linear3DFitterConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::Linear3DFitterConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLLinear3DFitter(void *p) {
      delete ((::larlitecv::Linear3DFitter*)p);
   }
   static void deleteArray_larlitecvcLcLLinear3DFitter(void *p) {
      delete [] ((::larlitecv::Linear3DFitter*)p);
   }
   static void destruct_larlitecvcLcLLinear3DFitter(void *p) {
      typedef ::larlitecv::Linear3DFitter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::Linear3DFitter

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStarDirNode(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarDirNode : new ::larlitecv::AStarDirNode;
   }
   static void *newArray_larlitecvcLcLAStarDirNode(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarDirNode[nElements] : new ::larlitecv::AStarDirNode[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarDirNode(void *p) {
      delete ((::larlitecv::AStarDirNode*)p);
   }
   static void deleteArray_larlitecvcLcLAStarDirNode(void *p) {
      delete [] ((::larlitecv::AStarDirNode*)p);
   }
   static void destruct_larlitecvcLcLAStarDirNode(void *p) {
      typedef ::larlitecv::AStarDirNode current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarDirNode

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStarDirAlgoConfig(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarDirAlgoConfig : new ::larlitecv::AStarDirAlgoConfig;
   }
   static void *newArray_larlitecvcLcLAStarDirAlgoConfig(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStarDirAlgoConfig[nElements] : new ::larlitecv::AStarDirAlgoConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarDirAlgoConfig(void *p) {
      delete ((::larlitecv::AStarDirAlgoConfig*)p);
   }
   static void deleteArray_larlitecvcLcLAStarDirAlgoConfig(void *p) {
      delete [] ((::larlitecv::AStarDirAlgoConfig*)p);
   }
   static void destruct_larlitecvcLcLAStarDirAlgoConfig(void *p) {
      typedef ::larlitecv::AStarDirAlgoConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarDirAlgoConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStarDirAlgo(void *p) {
      delete ((::larlitecv::AStarDirAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLAStarDirAlgo(void *p) {
      delete [] ((::larlitecv::AStarDirAlgo*)p);
   }
   static void destruct_larlitecvcLcLAStarDirAlgo(void *p) {
      typedef ::larlitecv::AStarDirAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStarDirAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStar3DNode(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStar3DNode : new ::larlitecv::AStar3DNode;
   }
   static void *newArray_larlitecvcLcLAStar3DNode(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStar3DNode[nElements] : new ::larlitecv::AStar3DNode[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStar3DNode(void *p) {
      delete ((::larlitecv::AStar3DNode*)p);
   }
   static void deleteArray_larlitecvcLcLAStar3DNode(void *p) {
      delete [] ((::larlitecv::AStar3DNode*)p);
   }
   static void destruct_larlitecvcLcLAStar3DNode(void *p) {
      typedef ::larlitecv::AStar3DNode current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStar3DNode

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLAStar3DAlgoConfig(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStar3DAlgoConfig : new ::larlitecv::AStar3DAlgoConfig;
   }
   static void *newArray_larlitecvcLcLAStar3DAlgoConfig(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::AStar3DAlgoConfig[nElements] : new ::larlitecv::AStar3DAlgoConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStar3DAlgoConfig(void *p) {
      delete ((::larlitecv::AStar3DAlgoConfig*)p);
   }
   static void deleteArray_larlitecvcLcLAStar3DAlgoConfig(void *p) {
      delete [] ((::larlitecv::AStar3DAlgoConfig*)p);
   }
   static void destruct_larlitecvcLcLAStar3DAlgoConfig(void *p) {
      typedef ::larlitecv::AStar3DAlgoConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStar3DAlgoConfig

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLAStar3DAlgo(void *p) {
      delete ((::larlitecv::AStar3DAlgo*)p);
   }
   static void deleteArray_larlitecvcLcLAStar3DAlgo(void *p) {
      delete [] ((::larlitecv::AStar3DAlgo*)p);
   }
   static void destruct_larlitecvcLcLAStar3DAlgo(void *p) {
      typedef ::larlitecv::AStar3DAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::AStar3DAlgo

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLBezierCurve(void *p) {
      delete ((::larlitecv::BezierCurve*)p);
   }
   static void deleteArray_larlitecvcLcLBezierCurve(void *p) {
      delete [] ((::larlitecv::BezierCurve*)p);
   }
   static void destruct_larlitecvcLcLBezierCurve(void *p) {
      typedef ::larlitecv::BezierCurve current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::BezierCurve

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
   static TClass *vectorlElarlitecvcLcLPointInfogR_Dictionary();
   static void vectorlElarlitecvcLcLPointInfogR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLPointInfogR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLPointInfogR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLPointInfogR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLPointInfogR(void *p);
   static void destruct_vectorlElarlitecvcLcLPointInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::PointInfo>*)
   {
      vector<larlitecv::PointInfo> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::PointInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::PointInfo>", -2, "vector", 210,
                  typeid(vector<larlitecv::PointInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLPointInfogR_Dictionary, isa_proxy, 0,
                  sizeof(vector<larlitecv::PointInfo>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLPointInfogR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLPointInfogR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLPointInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLPointInfogR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLPointInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::PointInfo> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::PointInfo>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLPointInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::PointInfo>*)0x0)->GetClass();
      vectorlElarlitecvcLcLPointInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLPointInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLPointInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::PointInfo> : new vector<larlitecv::PointInfo>;
   }
   static void *newArray_vectorlElarlitecvcLcLPointInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::PointInfo>[nElements] : new vector<larlitecv::PointInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLPointInfogR(void *p) {
      delete ((vector<larlitecv::PointInfo>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLPointInfogR(void *p) {
      delete [] ((vector<larlitecv::PointInfo>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLPointInfogR(void *p) {
      typedef vector<larlitecv::PointInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::PointInfo>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLBoundarySpacePointgR_Dictionary();
   static void vectorlElarlitecvcLcLBoundarySpacePointgR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLBoundarySpacePointgR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p);
   static void destruct_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::BoundarySpacePoint>*)
   {
      vector<larlitecv::BoundarySpacePoint> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::BoundarySpacePoint>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::BoundarySpacePoint>", -2, "vector", 210,
                  typeid(vector<larlitecv::BoundarySpacePoint>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLBoundarySpacePointgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<larlitecv::BoundarySpacePoint>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLBoundarySpacePointgR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLBoundarySpacePointgR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLBoundarySpacePointgR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLBoundarySpacePointgR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLBoundarySpacePointgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::BoundarySpacePoint> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::BoundarySpacePoint>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLBoundarySpacePointgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::BoundarySpacePoint>*)0x0)->GetClass();
      vectorlElarlitecvcLcLBoundarySpacePointgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLBoundarySpacePointgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BoundarySpacePoint> : new vector<larlitecv::BoundarySpacePoint>;
   }
   static void *newArray_vectorlElarlitecvcLcLBoundarySpacePointgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BoundarySpacePoint>[nElements] : new vector<larlitecv::BoundarySpacePoint>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p) {
      delete ((vector<larlitecv::BoundarySpacePoint>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p) {
      delete [] ((vector<larlitecv::BoundarySpacePoint>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLBoundarySpacePointgR(void *p) {
      typedef vector<larlitecv::BoundarySpacePoint> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::BoundarySpacePoint>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLBoundaryEndPtgR_Dictionary();
   static void vectorlElarlitecvcLcLBoundaryEndPtgR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLBoundaryEndPtgR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p);
   static void destruct_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::BoundaryEndPt>*)
   {
      vector<larlitecv::BoundaryEndPt> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::BoundaryEndPt>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::BoundaryEndPt>", -2, "vector", 210,
                  typeid(vector<larlitecv::BoundaryEndPt>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLBoundaryEndPtgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<larlitecv::BoundaryEndPt>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLBoundaryEndPtgR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLBoundaryEndPtgR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLBoundaryEndPtgR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLBoundaryEndPtgR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLBoundaryEndPtgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::BoundaryEndPt> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::BoundaryEndPt>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLBoundaryEndPtgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::BoundaryEndPt>*)0x0)->GetClass();
      vectorlElarlitecvcLcLBoundaryEndPtgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLBoundaryEndPtgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BoundaryEndPt> : new vector<larlitecv::BoundaryEndPt>;
   }
   static void *newArray_vectorlElarlitecvcLcLBoundaryEndPtgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BoundaryEndPt>[nElements] : new vector<larlitecv::BoundaryEndPt>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p) {
      delete ((vector<larlitecv::BoundaryEndPt>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p) {
      delete [] ((vector<larlitecv::BoundaryEndPt>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLBoundaryEndPtgR(void *p) {
      typedef vector<larlitecv::BoundaryEndPt> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::BoundaryEndPt>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLBMTrackCluster2DgR_Dictionary();
   static void vectorlElarlitecvcLcLBMTrackCluster2DgR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLBMTrackCluster2DgR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p);
   static void destruct_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::BMTrackCluster2D>*)
   {
      vector<larlitecv::BMTrackCluster2D> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::BMTrackCluster2D>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::BMTrackCluster2D>", -2, "vector", 210,
                  typeid(vector<larlitecv::BMTrackCluster2D>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLBMTrackCluster2DgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<larlitecv::BMTrackCluster2D>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLBMTrackCluster2DgR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLBMTrackCluster2DgR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLBMTrackCluster2DgR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLBMTrackCluster2DgR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLBMTrackCluster2DgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::BMTrackCluster2D> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::BMTrackCluster2D>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLBMTrackCluster2DgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::BMTrackCluster2D>*)0x0)->GetClass();
      vectorlElarlitecvcLcLBMTrackCluster2DgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLBMTrackCluster2DgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BMTrackCluster2D> : new vector<larlitecv::BMTrackCluster2D>;
   }
   static void *newArray_vectorlElarlitecvcLcLBMTrackCluster2DgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::BMTrackCluster2D>[nElements] : new vector<larlitecv::BMTrackCluster2D>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p) {
      delete ((vector<larlitecv::BMTrackCluster2D>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p) {
      delete [] ((vector<larlitecv::BMTrackCluster2D>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLBMTrackCluster2DgR(void *p) {
      typedef vector<larlitecv::BMTrackCluster2D> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::BMTrackCluster2D>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLAStarNodegR_Dictionary();
   static void vectorlElarlitecvcLcLAStarNodegR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLAStarNodegR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLAStarNodegR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLAStarNodegR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLAStarNodegR(void *p);
   static void destruct_vectorlElarlitecvcLcLAStarNodegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::AStarNode>*)
   {
      vector<larlitecv::AStarNode> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::AStarNode>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::AStarNode>", -2, "vector", 210,
                  typeid(vector<larlitecv::AStarNode>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLAStarNodegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<larlitecv::AStarNode>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLAStarNodegR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLAStarNodegR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLAStarNodegR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLAStarNodegR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLAStarNodegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::AStarNode> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::AStarNode>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLAStarNodegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::AStarNode>*)0x0)->GetClass();
      vectorlElarlitecvcLcLAStarNodegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLAStarNodegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLAStarNodegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStarNode> : new vector<larlitecv::AStarNode>;
   }
   static void *newArray_vectorlElarlitecvcLcLAStarNodegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStarNode>[nElements] : new vector<larlitecv::AStarNode>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLAStarNodegR(void *p) {
      delete ((vector<larlitecv::AStarNode>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLAStarNodegR(void *p) {
      delete [] ((vector<larlitecv::AStarNode>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLAStarNodegR(void *p) {
      typedef vector<larlitecv::AStarNode> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::AStarNode>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLAStarDirNodegR_Dictionary();
   static void vectorlElarlitecvcLcLAStarDirNodegR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLAStarDirNodegR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLAStarDirNodegR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLAStarDirNodegR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLAStarDirNodegR(void *p);
   static void destruct_vectorlElarlitecvcLcLAStarDirNodegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::AStarDirNode>*)
   {
      vector<larlitecv::AStarDirNode> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::AStarDirNode>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::AStarDirNode>", -2, "vector", 210,
                  typeid(vector<larlitecv::AStarDirNode>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLAStarDirNodegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<larlitecv::AStarDirNode>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLAStarDirNodegR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLAStarDirNodegR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLAStarDirNodegR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLAStarDirNodegR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLAStarDirNodegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::AStarDirNode> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::AStarDirNode>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLAStarDirNodegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::AStarDirNode>*)0x0)->GetClass();
      vectorlElarlitecvcLcLAStarDirNodegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLAStarDirNodegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLAStarDirNodegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStarDirNode> : new vector<larlitecv::AStarDirNode>;
   }
   static void *newArray_vectorlElarlitecvcLcLAStarDirNodegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStarDirNode>[nElements] : new vector<larlitecv::AStarDirNode>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLAStarDirNodegR(void *p) {
      delete ((vector<larlitecv::AStarDirNode>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLAStarDirNodegR(void *p) {
      delete [] ((vector<larlitecv::AStarDirNode>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLAStarDirNodegR(void *p) {
      typedef vector<larlitecv::AStarDirNode> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::AStarDirNode>

namespace ROOT {
   static TClass *vectorlElarlitecvcLcLAStar3DNodegR_Dictionary();
   static void vectorlElarlitecvcLcLAStar3DNodegR_TClassManip(TClass*);
   static void *new_vectorlElarlitecvcLcLAStar3DNodegR(void *p = 0);
   static void *newArray_vectorlElarlitecvcLcLAStar3DNodegR(Long_t size, void *p);
   static void delete_vectorlElarlitecvcLcLAStar3DNodegR(void *p);
   static void deleteArray_vectorlElarlitecvcLcLAStar3DNodegR(void *p);
   static void destruct_vectorlElarlitecvcLcLAStar3DNodegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlitecv::AStar3DNode>*)
   {
      vector<larlitecv::AStar3DNode> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlitecv::AStar3DNode>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlitecv::AStar3DNode>", -2, "vector", 210,
                  typeid(vector<larlitecv::AStar3DNode>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElarlitecvcLcLAStar3DNodegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<larlitecv::AStar3DNode>) );
      instance.SetNew(&new_vectorlElarlitecvcLcLAStar3DNodegR);
      instance.SetNewArray(&newArray_vectorlElarlitecvcLcLAStar3DNodegR);
      instance.SetDelete(&delete_vectorlElarlitecvcLcLAStar3DNodegR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecvcLcLAStar3DNodegR);
      instance.SetDestructor(&destruct_vectorlElarlitecvcLcLAStar3DNodegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlitecv::AStar3DNode> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlitecv::AStar3DNode>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecvcLcLAStar3DNodegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlitecv::AStar3DNode>*)0x0)->GetClass();
      vectorlElarlitecvcLcLAStar3DNodegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecvcLcLAStar3DNodegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecvcLcLAStar3DNodegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStar3DNode> : new vector<larlitecv::AStar3DNode>;
   }
   static void *newArray_vectorlElarlitecvcLcLAStar3DNodegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<larlitecv::AStar3DNode>[nElements] : new vector<larlitecv::AStar3DNode>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecvcLcLAStar3DNodegR(void *p) {
      delete ((vector<larlitecv::AStar3DNode>*)p);
   }
   static void deleteArray_vectorlElarlitecvcLcLAStar3DNodegR(void *p) {
      delete [] ((vector<larlitecv::AStar3DNode>*)p);
   }
   static void destruct_vectorlElarlitecvcLcLAStar3DNodegR(void *p) {
      typedef vector<larlitecv::AStar3DNode> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlitecv::AStar3DNode>

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
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = 0);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 518,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<bool>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bool>*)0x0)->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete ((vector<bool>*)p);
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] ((vector<bool>*)p);
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace ROOT {
   static TClass *setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_Dictionary();
   static void setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_TClassManip(TClass*);
   static void *new_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p = 0);
   static void *newArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(Long_t size, void *p);
   static void delete_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p);
   static void deleteArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p);
   static void destruct_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<larlitecv::AStarNode,larlitecv::node_pos_compare>*)
   {
      set<larlitecv::AStarNode,larlitecv::node_pos_compare> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<larlitecv::AStarNode,larlitecv::node_pos_compare>));
      static ::ROOT::TGenericClassInfo 
         instance("set<larlitecv::AStarNode,larlitecv::node_pos_compare>", -2, "set", 90,
                  typeid(set<larlitecv::AStarNode,larlitecv::node_pos_compare>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_Dictionary, isa_proxy, 0,
                  sizeof(set<larlitecv::AStarNode,larlitecv::node_pos_compare>) );
      instance.SetNew(&new_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR);
      instance.SetNewArray(&newArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR);
      instance.SetDelete(&delete_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR);
      instance.SetDeleteArray(&deleteArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR);
      instance.SetDestructor(&destruct_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<larlitecv::AStarNode,larlitecv::node_pos_compare> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const set<larlitecv::AStarNode,larlitecv::node_pos_compare>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<larlitecv::AStarNode,larlitecv::node_pos_compare>*)0x0)->GetClass();
      setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_TClassManip(theClass);
   return theClass;
   }

   static void setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<larlitecv::AStarNode,larlitecv::node_pos_compare> : new set<larlitecv::AStarNode,larlitecv::node_pos_compare>;
   }
   static void *newArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<larlitecv::AStarNode,larlitecv::node_pos_compare>[nElements] : new set<larlitecv::AStarNode,larlitecv::node_pos_compare>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p) {
      delete ((set<larlitecv::AStarNode,larlitecv::node_pos_compare>*)p);
   }
   static void deleteArray_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p) {
      delete [] ((set<larlitecv::AStarNode,larlitecv::node_pos_compare>*)p);
   }
   static void destruct_setlElarlitecvcLcLAStarNodecOlarlitecvcLcLnode_pos_comparegR(void *p) {
      typedef set<larlitecv::AStarNode,larlitecv::node_pos_compare> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<larlitecv::AStarNode,larlitecv::node_pos_compare>

namespace {
  void TriggerDictionaryInitialization_libLArliteCVThruMu_Impl() {
    static const char* headers[] = {
"AStar3DAlgo.h",
"AStar3DPostProcessor.h",
"AStarDirAlgo.h",
"AStarGridAlgo.h",
"BezierCurve.h",
"BMTrackCluster2D.h",
"BMTrackCluster3D.h",
"BoundaryEndPt.h",
"BoundaryIntersectionAlgo.h",
"BoundaryMatchAlgo.h",
"BoundaryMatchArrays.h",
"BoundaryMuonTaggerAlgo.h",
"BoundaryMuonTaggerTypes.h",
"BoundarySpacePoint.h",
"ConfigBoundaryMuonTaggerAlgo.h",
"DBClusterMergerAlgo.h",
"EndPointFilter.h",
"FlashMuonTaggerAlgo.h",
"FlashMuonTaggerConfig.h",
"Linear3DChargeTagger.h",
"Linear3DChargeTaggerTypes.h",
"Linear3DFitter.h",
"Linear3DPostProcessor_crbversion.h",
"Linear3DPostProcessor.h",
"LineRegionTest.h",
"TrackTests.h",
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
"/home/barnchri/CV_Directories_Test/larlite/core/../UserDev/BasicTool/GeoAlgo",
"/home/barnchri/CV_Directories_Test/larlitecv/build/include",
"/home/barnchri/CV_Directories_Test/LArCV/app/ann_1.1.2/include",
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/app/ThruMu/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVThruMu dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryMatchAlgo.h")))  BoundaryMatchArrays;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryMuonTaggerAlgo.h")))  ConfigBoundaryMuonTaggerAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DPostProcessor.h")))  BoundaryEndPt;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DPostProcessor.h")))  BoundarySpacePoint;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DPostProcessor.h")))  BMTrackCluster2D;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DPostProcessor.h")))  BMTrackCluster3D;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryIntersectionAlgo.h")))  BoundaryIntersectionAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryMuonTaggerAlgo.h")))  BoundaryMuonTaggerAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$FlashMuonTaggerAlgo.h")))  FlashMuonTaggerConfig;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$FlashMuonTaggerAlgo.h")))  FlashMuonTaggerAlgo;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarGridAlgo.h")))  AStarNode;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarGridAlgo.h")))  AStarSet;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarGridAlgo.h")))  AStarAlgoConfig;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarGridAlgo.h")))  AStarGridAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryMuonTaggerAlgo.h")))  PointInfo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BoundaryMuonTaggerAlgo.h")))  PointInfoList;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$Linear3DFitter.h")))  Linear3DFitterConfig;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$Linear3DFitter.h")))  Linear3DFitter;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarDirAlgo.h")))  AStarDirNode;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarDirAlgo.h")))  AStarDirAlgoConfig;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStarDirAlgo.h")))  AStarDirAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DAlgo.h")))  AStar3DNode;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DAlgo.h")))  AStar3DAlgoConfig;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$AStar3DAlgo.h")))  AStar3DAlgo;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$BezierCurve.h")))  BezierCurve;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVThruMu dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "AStar3DAlgo.h"
#include "AStar3DPostProcessor.h"
#include "AStarDirAlgo.h"
#include "AStarGridAlgo.h"
#include "BezierCurve.h"
#include "BMTrackCluster2D.h"
#include "BMTrackCluster3D.h"
#include "BoundaryEndPt.h"
#include "BoundaryIntersectionAlgo.h"
#include "BoundaryMatchAlgo.h"
#include "BoundaryMatchArrays.h"
#include "BoundaryMuonTaggerAlgo.h"
#include "BoundaryMuonTaggerTypes.h"
#include "BoundarySpacePoint.h"
#include "ConfigBoundaryMuonTaggerAlgo.h"
#include "DBClusterMergerAlgo.h"
#include "EndPointFilter.h"
#include "FlashMuonTaggerAlgo.h"
#include "FlashMuonTaggerConfig.h"
#include "Linear3DChargeTagger.h"
#include "Linear3DChargeTaggerTypes.h"
#include "Linear3DFitter.h"
#include "Linear3DPostProcessor_crbversion.h"
#include "Linear3DPostProcessor.h"
#include "LineRegionTest.h"
#include "TrackTests.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlitecv::AStar3DAlgo", payloadCode, "@",
"larlitecv::AStar3DAlgoConfig", payloadCode, "@",
"larlitecv::AStar3DNode", payloadCode, "@",
"larlitecv::AStarAlgoConfig", payloadCode, "@",
"larlitecv::AStarDirAlgo", payloadCode, "@",
"larlitecv::AStarDirAlgoConfig", payloadCode, "@",
"larlitecv::AStarDirNode", payloadCode, "@",
"larlitecv::AStarGridAlgo", payloadCode, "@",
"larlitecv::AStarNode", payloadCode, "@",
"larlitecv::AStarSet", payloadCode, "@",
"larlitecv::BMTrackCluster2D", payloadCode, "@",
"larlitecv::BMTrackCluster3D", payloadCode, "@",
"larlitecv::BezierCurve", payloadCode, "@",
"larlitecv::BoundaryEndPt", payloadCode, "@",
"larlitecv::BoundaryIntersectionAlgo", payloadCode, "@",
"larlitecv::BoundaryMatchArrays", payloadCode, "@",
"larlitecv::BoundaryMuonTaggerAlgo", payloadCode, "@",
"larlitecv::BoundarySpacePoint", payloadCode, "@",
"larlitecv::ConfigBoundaryMuonTaggerAlgo", payloadCode, "@",
"larlitecv::FlashMuonTaggerAlgo", payloadCode, "@",
"larlitecv::FlashMuonTaggerConfig", payloadCode, "@",
"larlitecv::Linear3DFitter", payloadCode, "@",
"larlitecv::Linear3DFitterConfig", payloadCode, "@",
"larlitecv::PointInfo", payloadCode, "@",
"larlitecv::PointInfoList", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVThruMu",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVThruMu_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVThruMu_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVThruMu() {
  TriggerDictionaryInitialization_libLArliteCVThruMu_Impl();
}
