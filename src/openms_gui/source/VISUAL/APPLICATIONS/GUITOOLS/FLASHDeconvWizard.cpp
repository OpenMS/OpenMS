// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

/**
  @page TOPP_FLASHDeconvWizard FLASHDeconvWizard

  @brief An assistant for FLASHDeconv execution.

  The implementation of FLASHDeconvWizard is heavily inspired by the SwathWizard.
  The Wizard helps the user to run FLASHDeconv for Top-down proteomics analysis. (@ref TOPP_FLASHDeconv tool)

  Users can enter the required input data (mzML MS/MS data) in dedicated fields, usually by drag'n'droping files from the
  operating systems' file explorer (Explorer, Nautilus, Finder...).
  The main output of the Wizard is deconvolved feature files (*.tsv) from FLASHDeconv. Optional output files are as follows:
    - deconvoluted MSn spectra files (*.tsv)
    - deconvoluted mzML spectra file (*.mzML)
    - deconvoluted MS1 in ProMex output format (*.ms1ft)
    - deconvoluted MSn spectra files in TopFD output format (*.msalign)
    - deconvoluted MS1 feature file in TopFD output format (*.feature)
*/

// QT
#include <QApplication>

// OpenMS
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/VISUAL/APPLICATIONS/FLASHDeconvWizardBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

// STL
#include <iostream>
#include <map>

#ifdef OPENMS_WINDOWSPLATFORM
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501 // Win XP (and above)
  #endif
  #include <Windows.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "FLASHDeconvWizard";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage(Logger::LogStream& stream = OpenMS_Log_info)
{
  stream << "\n"
         << tool_name << " -- An assistant for FLASHDeconv.\n"
         << "\n"
         << "Usage: \n"
         << " " << tool_name << " [options] [files]\n"
         << "\n"
         << "Options are:\n"
         << "  --help           Shows this help\n"
         << "  --debug          Enables debug messages\n"
         << endl;
}

int main(int argc, const char** argv)
{
  // list of all the valid options
  std::map<std::string, std::string> valid_options, valid_flags, option_lists;
  valid_flags["--help"] = "help";
  valid_flags["--debug"] = "debug";

  Param param;
  param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

  // '--help' given
  if (param.exists("help"))
  {
    print_usage();
    return 0;
  }

  // '-debug' given
  if (param.exists("debug"))
  {
    OPENMS_LOG_INFO << "Debug flag provided. Enabling 'OPENMS_LOG_DEBUG' ..." << std::endl;
    OpenMS_Log_debug.insert(cout); // allows to use OPENMS_LOG_DEBUG << "something" << std::endl;
  }

  // test if unknown options were given
  if (param.exists("unknown"))
  {
    // if packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if (!(String(param.getValue("unknown").toString()).hasSubstring("-psn") && !String(param.getValue("unknown").toString()).hasSubstring(", ")))
    {
      OPENMS_LOG_ERROR << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage(OpenMS_Log_error);
      return 1;
    }
  }

  QApplicationTOPP a(argc, const_cast<char**>(argv));
  a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

  FLASHDeconvWizardBase fw;
  fw.show();

#ifdef OPENMS_WINDOWSPLATFORM
  FreeConsole();     // get rid of console window at this point (we will not see any console output from this point on)
  AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

  int result = a.exec();

  return result;
}
