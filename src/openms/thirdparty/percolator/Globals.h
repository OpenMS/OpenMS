/*******************************************************************************
Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*******************************************************************************/
// Static stub replacing Percolator's CMake-generated Globals.h (Globals.h.cmake
// in upstream, with the values of upstream's CommonCMake.txt substituted).
#ifndef GLOBALS_H_
#define GLOBALS_H_

#ifndef PIN_VERSION_MAJOR
  #define PIN_VERSION_MAJOR "1"
#endif

#ifndef POUT_VERSION_MAJOR
  #define POUT_VERSION_MAJOR "1"
#endif

#ifndef PIN_VERSION_MINOR
  #define PIN_VERSION_MINOR "3"
#endif

#ifndef POUT_VERSION_MINOR
  #define POUT_VERSION_MINOR "5"
#endif

#ifndef M_PI // for some reason Windows removed the Pi constant from math.h
  #define M_PI 3.14159265358979323846
#endif

#ifdef WIN32
  #include <float.h>
#endif

#ifdef _WIN32
  #ifndef __stdcall
    //#define __stdcall
  #endif
#endif
#ifdef _WIN32
  #include <fcntl.h>
#endif

#ifdef _MSC_VER
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
#endif



#ifdef __APPLE__
  #ifdef __INTEL_COMPILER
    #define isnan(x) _isnan(x)
    #define isinf(x) _isinf(x)
//  #else
//    #define isinf(x) std::isinf(x)
//    #define isnan(x) x!=x
  #endif
#endif


#define VERB (Globals::getInstance()->getVerbose())
#define NO_TERMINATE (Globals::getInstance()->getNoTerminate())
#include <time.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "MyException.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;

/* Process-wide settings: verbosity and the no-terminate flag. */
class Globals {
  public:
    static Globals* getInstance();
    // restores the default settings
    static void clean();
    
    int getVerbose() {
      return verbose;
    }
    void setVerbose(int verb) {
      verbose = verb;
    }
    void decVerbose() {
      verbose--;
    }
    void incVerbose() {
      verbose++;
    }
    void setNoTerminate(bool noTerminate) {
      noTerminate_ = noTerminate;
    }
    bool getNoTerminate() {
      return noTerminate_;
    }
    // For an error that --no-terminate lets the run survive. Without the
    // flag, throws MyException(error + thrownSuffix). With it, prints error
    // and recovery ("No-terminate flag set: ...") to stderr and returns; the
    // caller then does what the recovery text says.
    static void failUnlessNoTerminate(const std::string& error,
                                      const char* recovery,
                                      const char* thrownSuffix = "");
  private:
    Globals() = default;
    int verbose = 2;
    bool noTerminate_ = false;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*GLOBALS_H_*/

