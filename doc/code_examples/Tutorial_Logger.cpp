// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

// just a placeholder for the computations you would do in real code
void someFunction()
{
}

int main(int argc, const char** argv)
{
  //! [doxygen_snippet_Logger]

  ProgressLogger progresslogger;
  progresslogger.setLogType(ProgressLogger::CMD); // output to the terminal (std::cout)
  // Note: within a TOPP tool, you can use
  // progresslogger.setLogType(TOPPBase::log_type_);
  // to set the log-type (automatically set via commandline options)

  const int progress_steps = 200;
  // set start progress (0) and end (ms_run.size() = the number of spectra)
  progresslogger.startProgress(0, progress_steps, "Doing some calculation...");

  for (int i = 0; i < progress_steps; ++i) // in real code, iterate over some datastructure, e.g. an MSExperiments' spectra
  {
    // update progress
    progresslogger.setProgress(i);
    // do the actual calculations and processing ...
    someFunction();
  }

  progresslogger.endProgress();
//! [doxygen_snippet_Logger]

}
