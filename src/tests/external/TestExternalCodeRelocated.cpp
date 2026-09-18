// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/VersionInfo.h>

#include <iostream>
#include <string>

// Consumer program of an OpenMS installation whose prefix was renamed after
// `cmake --install` (relocated/CMakeLists.txt). Deliberately free of the OpenMS
// shared data: the installation this is built against leaves out the share
// component, so nothing here may reach for <prefix>/share/OpenMS.

using namespace OpenMS;

int main()
{
  // Compiled against the headers of the moved prefix and linked against the
  // libOpenMS installed there. Reaching main() at all means the loader resolved
  // that library from where OpenMSConfig.cmake said it was, i.e. the relocated
  // installation carries an RPATH (an entry on PATH on Windows) that survived
  // the move.
  const std::string version = VersionInfo::getVersion();
  std::cout << "Loaded OpenMS " << version << " from the relocated installation" << std::endl;

  // The library that was loaded has to be the one the package describes, not
  // another OpenMS that happened to be earlier on the search path. Compared
  // through the parsed struct rather than as a string: a pre-release build
  // appends "-pre-<identifier>" to the library's version, while the package
  // reports the plain major.minor.patch in OPENMS_EXPECTED_VERSION.
  const VersionInfo::VersionDetails expected = VersionInfo::VersionDetails::create(OPENMS_EXPECTED_VERSION);
  if (expected == VersionInfo::VersionDetails::EMPTY)
  {
    std::cerr << "cannot parse the version the package reported: " << OPENMS_EXPECTED_VERSION << "\n";
    return 1;
  }
  const VersionInfo::VersionDetails loaded = VersionInfo::getVersionStruct();
  if (loaded.version_major != expected.version_major
      || loaded.version_minor != expected.version_minor
      || loaded.version_patch != expected.version_patch)
  {
    std::cerr << "the package reports OpenMS " << OPENMS_EXPECTED_VERSION
              << ", but the library loaded from it is " << version << "\n";
    return 1;
  }

  std::cout << "All good and well!" << std::endl;
  return 0;
}
