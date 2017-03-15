// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIBasedIBaseDict

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
#include "DataCoordinator.h"
#include "FileManager.h"
#include "LarcvFileManager.h"
#include "LarliteFileManager.h"

// Header files passed via #pragma extra_include

namespace larlitecv {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *larlitecv_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("larlitecv", 0 /*version*/, "DataCoordinator.h", 20,
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
   static TClass *larlitecvcLcLLarliteFileManager_Dictionary();
   static void larlitecvcLcLLarliteFileManager_TClassManip(TClass*);
   static void delete_larlitecvcLcLLarliteFileManager(void *p);
   static void deleteArray_larlitecvcLcLLarliteFileManager(void *p);
   static void destruct_larlitecvcLcLLarliteFileManager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::LarliteFileManager*)
   {
      ::larlitecv::LarliteFileManager *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::LarliteFileManager));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::LarliteFileManager", "LarliteFileManager.h", 11,
                  typeid(::larlitecv::LarliteFileManager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLLarliteFileManager_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::LarliteFileManager) );
      instance.SetDelete(&delete_larlitecvcLcLLarliteFileManager);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLLarliteFileManager);
      instance.SetDestructor(&destruct_larlitecvcLcLLarliteFileManager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::LarliteFileManager*)
   {
      return GenerateInitInstanceLocal((::larlitecv::LarliteFileManager*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::LarliteFileManager*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLLarliteFileManager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::LarliteFileManager*)0x0)->GetClass();
      larlitecvcLcLLarliteFileManager_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLLarliteFileManager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLLarcvFileManager_Dictionary();
   static void larlitecvcLcLLarcvFileManager_TClassManip(TClass*);
   static void delete_larlitecvcLcLLarcvFileManager(void *p);
   static void deleteArray_larlitecvcLcLLarcvFileManager(void *p);
   static void destruct_larlitecvcLcLLarcvFileManager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::LarcvFileManager*)
   {
      ::larlitecv::LarcvFileManager *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::LarcvFileManager));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::LarcvFileManager", "LarcvFileManager.h", 11,
                  typeid(::larlitecv::LarcvFileManager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLLarcvFileManager_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::LarcvFileManager) );
      instance.SetDelete(&delete_larlitecvcLcLLarcvFileManager);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLLarcvFileManager);
      instance.SetDestructor(&destruct_larlitecvcLcLLarcvFileManager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::LarcvFileManager*)
   {
      return GenerateInitInstanceLocal((::larlitecv::LarcvFileManager*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::LarcvFileManager*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLLarcvFileManager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::LarcvFileManager*)0x0)->GetClass();
      larlitecvcLcLLarcvFileManager_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLLarcvFileManager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecvcLcLDataCoordinator_Dictionary();
   static void larlitecvcLcLDataCoordinator_TClassManip(TClass*);
   static void *new_larlitecvcLcLDataCoordinator(void *p = 0);
   static void *newArray_larlitecvcLcLDataCoordinator(Long_t size, void *p);
   static void delete_larlitecvcLcLDataCoordinator(void *p);
   static void deleteArray_larlitecvcLcLDataCoordinator(void *p);
   static void destruct_larlitecvcLcLDataCoordinator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlitecv::DataCoordinator*)
   {
      ::larlitecv::DataCoordinator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlitecv::DataCoordinator));
      static ::ROOT::TGenericClassInfo 
         instance("larlitecv::DataCoordinator", "DataCoordinator.h", 24,
                  typeid(::larlitecv::DataCoordinator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecvcLcLDataCoordinator_Dictionary, isa_proxy, 4,
                  sizeof(::larlitecv::DataCoordinator) );
      instance.SetNew(&new_larlitecvcLcLDataCoordinator);
      instance.SetNewArray(&newArray_larlitecvcLcLDataCoordinator);
      instance.SetDelete(&delete_larlitecvcLcLDataCoordinator);
      instance.SetDeleteArray(&deleteArray_larlitecvcLcLDataCoordinator);
      instance.SetDestructor(&destruct_larlitecvcLcLDataCoordinator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlitecv::DataCoordinator*)
   {
      return GenerateInitInstanceLocal((::larlitecv::DataCoordinator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlitecv::DataCoordinator*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecvcLcLDataCoordinator_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlitecv::DataCoordinator*)0x0)->GetClass();
      larlitecvcLcLDataCoordinator_TClassManip(theClass);
   return theClass;
   }

   static void larlitecvcLcLDataCoordinator_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLLarliteFileManager(void *p) {
      delete ((::larlitecv::LarliteFileManager*)p);
   }
   static void deleteArray_larlitecvcLcLLarliteFileManager(void *p) {
      delete [] ((::larlitecv::LarliteFileManager*)p);
   }
   static void destruct_larlitecvcLcLLarliteFileManager(void *p) {
      typedef ::larlitecv::LarliteFileManager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::LarliteFileManager

namespace ROOT {
   // Wrapper around operator delete
   static void delete_larlitecvcLcLLarcvFileManager(void *p) {
      delete ((::larlitecv::LarcvFileManager*)p);
   }
   static void deleteArray_larlitecvcLcLLarcvFileManager(void *p) {
      delete [] ((::larlitecv::LarcvFileManager*)p);
   }
   static void destruct_larlitecvcLcLLarcvFileManager(void *p) {
      typedef ::larlitecv::LarcvFileManager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::LarcvFileManager

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecvcLcLDataCoordinator(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::DataCoordinator : new ::larlitecv::DataCoordinator;
   }
   static void *newArray_larlitecvcLcLDataCoordinator(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlitecv::DataCoordinator[nElements] : new ::larlitecv::DataCoordinator[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecvcLcLDataCoordinator(void *p) {
      delete ((::larlitecv::DataCoordinator*)p);
   }
   static void deleteArray_larlitecvcLcLDataCoordinator(void *p) {
      delete [] ((::larlitecv::DataCoordinator*)p);
   }
   static void destruct_larlitecvcLcLDataCoordinator(void *p) {
      typedef ::larlitecv::DataCoordinator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlitecv::DataCoordinator

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 210,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *maplEstringcOvectorlEstringgRsPgR_Dictionary();
   static void maplEstringcOvectorlEstringgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlEstringgRsPgR(void *p = 0);
   static void *newArray_maplEstringcOvectorlEstringgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlEstringgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlEstringgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlEstringgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<string> >*)
   {
      map<string,vector<string> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<string> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<string> >", -2, "map", 96,
                  typeid(map<string,vector<string> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlEstringgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,vector<string> >) );
      instance.SetNew(&new_maplEstringcOvectorlEstringgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlEstringgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlEstringgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<string> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<string,vector<string> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlEstringgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,vector<string> >*)0x0)->GetClass();
      maplEstringcOvectorlEstringgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlEstringgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlEstringgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<string> > : new map<string,vector<string> >;
   }
   static void *newArray_maplEstringcOvectorlEstringgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<string> >[nElements] : new map<string,vector<string> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlEstringgRsPgR(void *p) {
      delete ((map<string,vector<string> >*)p);
   }
   static void deleteArray_maplEstringcOvectorlEstringgRsPgR(void *p) {
      delete [] ((map<string,vector<string> >*)p);
   }
   static void destruct_maplEstringcOvectorlEstringgRsPgR(void *p) {
      typedef map<string,vector<string> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<string> >

namespace ROOT {
   static TClass *maplEstringcOstringgR_Dictionary();
   static void maplEstringcOstringgR_TClassManip(TClass*);
   static void *new_maplEstringcOstringgR(void *p = 0);
   static void *newArray_maplEstringcOstringgR(Long_t size, void *p);
   static void delete_maplEstringcOstringgR(void *p);
   static void deleteArray_maplEstringcOstringgR(void *p);
   static void destruct_maplEstringcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,string>*)
   {
      map<string,string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,string>", -2, "map", 96,
                  typeid(map<string,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOstringgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,string>) );
      instance.SetNew(&new_maplEstringcOstringgR);
      instance.SetNewArray(&newArray_maplEstringcOstringgR);
      instance.SetDelete(&delete_maplEstringcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOstringgR);
      instance.SetDestructor(&destruct_maplEstringcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<string,string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,string>*)0x0)->GetClass();
      maplEstringcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,string> : new map<string,string>;
   }
   static void *newArray_maplEstringcOstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,string>[nElements] : new map<string,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOstringgR(void *p) {
      delete ((map<string,string>*)p);
   }
   static void deleteArray_maplEstringcOstringgR(void *p) {
      delete [] ((map<string,string>*)p);
   }
   static void destruct_maplEstringcOstringgR(void *p) {
      typedef map<string,string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,string>

namespace ROOT {
   static TClass *maplEstringcOlarlitecvcLcLFileManagermUgR_Dictionary();
   static void maplEstringcOlarlitecvcLcLFileManagermUgR_TClassManip(TClass*);
   static void *new_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p = 0);
   static void *newArray_maplEstringcOlarlitecvcLcLFileManagermUgR(Long_t size, void *p);
   static void delete_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p);
   static void deleteArray_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p);
   static void destruct_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,larlitecv::FileManager*>*)
   {
      map<string,larlitecv::FileManager*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,larlitecv::FileManager*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,larlitecv::FileManager*>", -2, "map", 96,
                  typeid(map<string,larlitecv::FileManager*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOlarlitecvcLcLFileManagermUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,larlitecv::FileManager*>) );
      instance.SetNew(&new_maplEstringcOlarlitecvcLcLFileManagermUgR);
      instance.SetNewArray(&newArray_maplEstringcOlarlitecvcLcLFileManagermUgR);
      instance.SetDelete(&delete_maplEstringcOlarlitecvcLcLFileManagermUgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOlarlitecvcLcLFileManagermUgR);
      instance.SetDestructor(&destruct_maplEstringcOlarlitecvcLcLFileManagermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,larlitecv::FileManager*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<string,larlitecv::FileManager*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOlarlitecvcLcLFileManagermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,larlitecv::FileManager*>*)0x0)->GetClass();
      maplEstringcOlarlitecvcLcLFileManagermUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOlarlitecvcLcLFileManagermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,larlitecv::FileManager*> : new map<string,larlitecv::FileManager*>;
   }
   static void *newArray_maplEstringcOlarlitecvcLcLFileManagermUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,larlitecv::FileManager*>[nElements] : new map<string,larlitecv::FileManager*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p) {
      delete ((map<string,larlitecv::FileManager*>*)p);
   }
   static void deleteArray_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p) {
      delete [] ((map<string,larlitecv::FileManager*>*)p);
   }
   static void destruct_maplEstringcOlarlitecvcLcLFileManagermUgR(void *p) {
      typedef map<string,larlitecv::FileManager*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,larlitecv::FileManager*>

namespace ROOT {
   static TClass *maplEstringcOintgR_Dictionary();
   static void maplEstringcOintgR_TClassManip(TClass*);
   static void *new_maplEstringcOintgR(void *p = 0);
   static void *newArray_maplEstringcOintgR(Long_t size, void *p);
   static void delete_maplEstringcOintgR(void *p);
   static void deleteArray_maplEstringcOintgR(void *p);
   static void destruct_maplEstringcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,int>*)
   {
      map<string,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,int>", -2, "map", 96,
                  typeid(map<string,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,int>) );
      instance.SetNew(&new_maplEstringcOintgR);
      instance.SetNewArray(&newArray_maplEstringcOintgR);
      instance.SetDelete(&delete_maplEstringcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOintgR);
      instance.SetDestructor(&destruct_maplEstringcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<string,int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,int>*)0x0)->GetClass();
      maplEstringcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,int> : new map<string,int>;
   }
   static void *newArray_maplEstringcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,int>[nElements] : new map<string,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOintgR(void *p) {
      delete ((map<string,int>*)p);
   }
   static void deleteArray_maplEstringcOintgR(void *p) {
      delete [] ((map<string,int>*)p);
   }
   static void destruct_maplEstringcOintgR(void *p) {
      typedef map<string,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,int>

namespace {
  void TriggerDictionaryInitialization_libLArliteCVBase_Impl() {
    static const char* headers[] = {
"DataCoordinator.h",
"FileManager.h",
"LarcvFileManager.h",
"LarliteFileManager.h",
0
    };
    static const char* includePaths[] = {
"/home/barnchri/CV_Directories_Test/larlite/core",
"/home/barnchri/CV_Directories_Test/LArCV/build/include",
"/usr/local/include",
"/usr/include/python2.7",
"/usr/include/x86_64-linux-gnu/python2.7",
"/usr/local/lib/python2.7/dist-packages/numpy/core/include",
"/home/barnchri/CV_Directories_Test/larlitecv/build/include",
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/core/Base/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVBase dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$LarliteFileManager.h")))  LarliteFileManager;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$LarcvFileManager.h")))  LarcvFileManager;}
namespace larlitecv{class __attribute__((annotate("$clingAutoload$DataCoordinator.h")))  DataCoordinator;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVBase dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "DataCoordinator.h"
#include "FileManager.h"
#include "LarcvFileManager.h"
#include "LarliteFileManager.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlitecv::DataCoordinator", payloadCode, "@",
"larlitecv::LarcvFileManager", payloadCode, "@",
"larlitecv::LarliteFileManager", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVBase",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVBase_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVBase_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVBase() {
  TriggerDictionaryInitialization_libLArliteCVBase_Impl();
}
