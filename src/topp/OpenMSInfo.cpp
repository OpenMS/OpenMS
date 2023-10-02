// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/BuildInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#ifdef _OPENMP
  #include "omp.h"
#endif

#include <iostream>

using namespace OpenMS;
using namespace std;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int main(int /*argc*/, const char ** /*argv*/)
{
  cout << "OpenMS Version:" << "\n";
  cout << "==================" << "\n";
  cout << "Version      : " << VersionInfo::getVersion() << "\n";
  cout << "Build time   : " << VersionInfo::getTime() << "\n";
  cout << "Git sha1     : " << VersionInfo::getRevision() << "\n";
  cout << "Git branch   : " << VersionInfo::getBranch() << "\n";
  cout << "\n";
  cout << "Installation information:" << "\n";
  cout << "==================" << "\n";
  cout << "Data path    : " << File::getOpenMSDataPath() << "\n";
  cout << "Temp path    : " << File::getTempDirectory() << "\n";
  cout << "Userdata path: " << File::getUserDirectory() << "\n";

  cout << "\n";
  cout << "Build information:" << "\n";
  cout << "==================" << "\n";
  cout << "Source path  : " << OPENMS_SOURCE_PATH << "\n";
  cout << "Binary path  : " << OPENMS_BINARY_PATH << "\n";
  cout << "Binary arch  : " << Internal::OpenMSOSInfo::getBinaryArchitecture() << "\n";
  cout << "Build type   : " << Internal::OpenMSBuildInfo::getBuildType() << "\n";
  #ifdef _OPENMP
  cout << "OpenMP       : " << "enabled (maxThreads = " << Internal::OpenMSBuildInfo::getOpenMPMaxNumThreads() << ")" << "\n";
  #else
  cout << "OpenMP       : " << "disabled" << "\n";
  #endif
  cout << "SIMD extensions : " << Internal::OpenMSOSInfo::getActiveSIMDExtensions() << "\n";
  cout << "\n";

  Internal::OpenMSOSInfo info = Internal::OpenMSOSInfo::getOSInfo();

  cout << "OS Information:" << "\n";
  cout << "==================" << "\n";
  cout << "Name: " << info.getOSAsString() << "\n";
  cout << "Version: " << info.getOSVersionAsString() << "\n";
  cout << "Architecture: " << info.getArchAsString() << "\n";
  cout << "\n";


  return 0;
}

/// @endcond
