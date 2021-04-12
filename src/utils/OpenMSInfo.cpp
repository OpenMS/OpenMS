// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
