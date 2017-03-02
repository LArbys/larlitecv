// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbarnchridICV_Directories_TestdIlarlitecvdIbuilddIHashlib2plusdIHashlib2plusDict

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
#include "hashlibpp.h"
#include "hl_exception.h"
#include "hl_hashwrapper.h"
#include "hl_md5.h"
#include "hl_md5wrapper.h"
#include "hl_sha1.h"
#include "hl_sha1wrapper.h"
#include "hl_sha256.h"
#include "hl_sha256wrapper.h"
#include "hl_sha2ext.h"
#include "hl_sha2mac.h"
#include "hl_sha384wrapper.h"
#include "hl_sha512wrapper.h"
#include "hl_types.h"
#include "hl_wrapperfactory.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_libLArliteCVHashlib2plus_Impl() {
    static const char* headers[] = {
"hashlibpp.h",
"hl_exception.h",
"hl_hashwrapper.h",
"hl_md5.h",
"hl_md5wrapper.h",
"hl_sha1.h",
"hl_sha1wrapper.h",
"hl_sha256.h",
"hl_sha256wrapper.h",
"hl_sha2ext.h",
"hl_sha2mac.h",
"hl_sha384wrapper.h",
"hl_sha512wrapper.h",
"hl_types.h",
"hl_wrapperfactory.h",
0
    };
    static const char* includePaths[] = {
"/home/spitzj/root-6.06.04/include",
"/home/barnchri/CV_Directories_Test/larlitecv/core/Hashlib2plus/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLArliteCVHashlib2plus dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLArliteCVHashlib2plus dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "hashlibpp.h"
#include "hl_exception.h"
#include "hl_hashwrapper.h"
#include "hl_md5.h"
#include "hl_md5wrapper.h"
#include "hl_sha1.h"
#include "hl_sha1wrapper.h"
#include "hl_sha256.h"
#include "hl_sha256wrapper.h"
#include "hl_sha2ext.h"
#include "hl_sha2mac.h"
#include "hl_sha384wrapper.h"
#include "hl_sha512wrapper.h"
#include "hl_types.h"
#include "hl_wrapperfactory.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLArliteCVHashlib2plus",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLArliteCVHashlib2plus_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLArliteCVHashlib2plus_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLArliteCVHashlib2plus() {
  TriggerDictionaryInitialization_libLArliteCVHashlib2plus_Impl();
}
