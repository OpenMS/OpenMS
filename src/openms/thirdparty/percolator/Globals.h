/*******************************************************************************
Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

*******************************************************************************/
// Static stub replacing Percolator's CMake-generated Globals.h (was Globals.h.cmake
// in upstream). Placeholder values substituted for paths we don't exercise in
// the vendored PSM-only subset.
#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <time.h>
#include <string>
#include <iosfwd>
#include <fstream>

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
#ifndef PERCOLATOR_IN_NAMESPACE
  #define PERCOLATOR_IN_NAMESPACE ""
#endif
#ifndef PERCOLATOR_OUT_NAMESPACE
  #define PERCOLATOR_OUT_NAMESPACE ""
#endif
#ifndef WRITABLE_DIR
  #define WRITABLE_DIR ""
#endif
#ifndef PIN_SCHEMA_LOCATION
  #define PIN_SCHEMA_LOCATION ""
#endif
#ifndef POUT_SCHEMA_LOCATION
  #define POUT_SCHEMA_LOCATION ""
#endif
#ifndef ELUDE_MODELS_PATH
  #define ELUDE_MODELS_PATH ""
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#ifdef WIN32
  #define C_DARRAY(name,nelem) double *name = (double *) _malloca((nelem) * sizeof(double));
  #define D_DARRAY(name) _freea(name);
  #include <float.h>
#else
  #define C_DARRAY(name,nelem) double name[nelem];
  #define D_DARRAY(name)
#endif

#ifdef _MSC_VER
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  // Modern MSVC has std::isfinite from <cmath>; the legacy
  // `#define isfinite _finite` macro that used to live here
  // pollutes <format> and Eigen (both reference std::isfinite /
  // std::_finite via their own helpers, and the function-like
  // macro mangles the std-qualified usage). Use std::isfinite at
  // call sites instead.
#endif

#include "Logger.h"
#include "MyException.h"

namespace OpenMS { namespace Internal { namespace Percolator {
using namespace std;

#define VERB (Globals::getInstance()->getVerbose())
#define NO_TERMINATE (Globals::getInstance()->getNoTerminate())
#define STD_CERR (*(Globals::getInstance()->getLogger()))

class Globals {
  public:
    virtual ~Globals();
    static Globals* getInstance();
    static void clean();
    Logger* getLogger();

    const std::string& getLogFile(){ return fileLog; }
    int getVerbose() { return verbose; }
    void setVerbose(int verb) { verbose = verb; }
    void decVerbose() { verbose--; }
    void incVerbose() { verbose++; }
    void setNoTerminate(bool noTerminate) { noTerminate_ = noTerminate; }
    bool getNoTerminate() { return noTerminate_; }
    void setLogFile(const std::string& filename);
    void initLogger();
    int redirectBuffer();
    void unredirectBuffer();
  private:
    Globals();
    int verbose;
    bool noTerminate_;
    static Globals* glob;
    Logger *log;
    std::string fileLog;
    bool buffer_redirected;
    std::streambuf* save_sbuf_cerr;
    std::ofstream ferr;
};

}}}  // namespace OpenMS::Internal::Percolator
#endif /*GLOBALS_H_*/
