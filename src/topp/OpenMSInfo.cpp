// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Colorizer.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/IndentedStream.h>
#include <OpenMS/SYSTEM/BuildInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/config.h>
#include <OpenMS/openms_data_path.h>

#ifdef _OPENMP
  #include "omp.h"
#endif

#include <iostream>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_OpenMSInfo OpenMSInfo

@brief Prints configurations details of %OpenMS (Version, Git hash, SIMD extensions, Multithreading), along with directories where auxilliary data
like modifications (UniMOD), Enzymes etc are taken from.

Some path's can be manipulated by the user by setting environment variables. If not set, the values are taken from the system defaults.

<ul>
  <li> <b>Data path:</b> controlled by the environment variable 'OPENMS_DATA_PATH'; the value should point to a %OpenMS share directory, e.g.
'c:/program files/OpenMS3.1/share/OpenMS' <li> <b>Temp path:</b> controlled by the environment variable 'OPENMS_TMPDIR'; the value should point to
where you want %OpenMS to store temporary data. <li> <b>Userdata path:</b> controlled by the environment variable 'OPENMS_HOME_PATH'; the value should
point to where you want %OpenMS to store user-realted data, e.g. the .OpenMS.ini.
</ul>

<B>This tool does not need/use any command line parameters.</B>

Example output:
@code

Full documentation: http://www.openms.de/doxygen/nightly/html/TOPP_OpenMSInfo.html
To cite OpenMS:
 + Pfeuffer, J., Bielow, C., Wein, S. et al.. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry
   data. Nat Methods (2024). doi:10.1038/s41592-024-02197-7.

<< OpenMS Version >>
Version      : 3.2.0
Build time   : Sep  18 2024, 14:14:53
Git sha1     : disabled
Git branch   : disabled

<< Installation information >>
Data path    : C:/dev/openms2/share/OpenMS
Temp path    : C:/Users/bielow/AppData/Local/Temp
Userdata path: C:/Users/bielow/

<< Build information >>
Source path  : C:/dev/openms/src/openms
Binary path  : C:/dev/openms_build/src/openms
Binary arch  : 64 bit
Build type   : Debug
LP-Solver    : COIN-OR
OpenMP       : enabled (maxThreads = 32)
SIMD extensions : SSE, SSE2, SSE3, SSE4.1, SSE4.2, AVX

<< OS Information >>
Name: Windows
Version: 10
Architecture: 64 bit

@endcode

*/
// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

/// this needs to be done before TOPPBase is initialized, since it will set OMP's max_threads to 1
const auto max_threads = Internal::OpenMSBuildInfo::getOpenMPMaxNumThreads();

class TOPPOpenMSInfo : public TOPPBase
{
public:
  TOPPOpenMSInfo(): TOPPBase("OpenMSInfo", "Prints configurations details of OpenMS.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerFlag_("p", "Print information (flag can also be omitted)");
    registerInputFile_("dummy", "<ignored>", "", "A fake input file, which is needed for some workflow systems to call this tool", false, true);
  }

  // Param getSubsectionDefaults_(const String& /*section*/) const override
  //{
  //   return SpectraMerger().getParameters();
  // }

  ExitCodes main_(int, const char**) override
  {
    IndentedStream is(cout, 0, 10);
    is << '\n'
       << bright("Full documentation: ") // the space is needed, otherwise the remaining line will be underlined on Windows..
       << underline(TOPPBase::getDocumentationURL()) << " " // the space is needed ...
       << '\n'
       << bright("To cite OpenMS:\n") << " + "
       << is.indent(3) << cite_openms.toString() << is.indent(0);

    is << "\n\n"
       << green("<< OpenMS Version >>\n")
       << "Version      : " << VersionInfo::getVersion() << "\n"
       << "Build time   : " << VersionInfo::getTime() << "\n"
       << "Git sha1     : " << VersionInfo::getRevision() << "\n"
       << "Git branch   : " << VersionInfo::getBranch() << "\n"
       << "\n"
       << green("<< Installation information >>\n")
       << "Data path    : " << File::getOpenMSDataPath() << "\n"
       << "Temp path    : " << File::getTempDirectory() << "\n"
       << "Userdata path: " << File::getUserDirectory() << "\n"
       << "\n"
       << green("<< Build information >>\n") 
       << "Source path  : " << OPENMS_SOURCE_PATH << "\n"
       << "Binary path  : " << OPENMS_BINARY_PATH << "\n"
       << "Binary arch  : " << Internal::OpenMSOSInfo::getBinaryArchitecture() << "\n"
       << "Build type   : " << Internal::OpenMSBuildInfo::getBuildType() << "\n";
#ifdef OPENMS_HAS_COINOR
    is << "LP-Solver    : COIN-OR\n";
#else
    cout << "LP-Solver    : GLPK\n";
#endif
#ifdef _OPENMP
    is << "OpenMP       : "
       << "enabled (maxThreads = " << max_threads << ")"
       << "\n";
#else
    is << "OpenMP       : "
       << "disabled"
       << "\n";
#endif
    is << "SIMD extensions : " << Internal::OpenMSOSInfo::getActiveSIMDExtensions() << "\n"
       << "\n";

    Internal::OpenMSOSInfo info = Internal::OpenMSOSInfo::getOSInfo();

    is << green("<< OS Information >>\n")
       << "Name: " << info.getOSAsString() << "\n"
       << "Version: " << info.getOSVersionAsString() << "\n"
       << "Architecture: " << info.getArchAsString() << "\n"
       << "\n";

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenMSInfo tool;
  // TOPPBase automatically shows the Help page if a TOPP tool is called without any parameters.
  // This tool is special: we want to print stuff in this case. So we pass a "-p" flag.
  const char* override_params[] = {argv[0], "-p"};
  if (argc == 1) // no args given
  {
    return tool.main(2, override_params);
  }
  return tool.main(argc, argv);
}

/// @endcond
