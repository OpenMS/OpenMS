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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QFile>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page

    @brief

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MzMLToMatlab.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MzMLToMatlab.html
*/

/// @cond TOPPCLASSES

class TOPPFLASHDeconvSpectrumSimulator :
    public TOPPBase
{
public:
  TOPPFLASHDeconvSpectrumSimulator() :
      TOPPBase("FLASHDeconvSpectrumSimulator", "Generate simulated top-down spectrum", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerIntOption_("in", "<number>", 10, "Input number of masses");
    registerDoubleOption_("min_mass", "<number>", 500, "min mass");
    registerDoubleOption_("max_mass", "<number>", 1e5, "max mass");

    registerDoubleOption_("min_mz", "<number>", 400, "min mz");
    registerDoubleOption_("max_mz", "<number>", 2000, "max mz");

    registerDoubleOption_("snr", "<number>", .1, "SNR");
    registerOutputFile_("out", "<mzML file>", "", "out mzML file");
    //setValidFormats_("in", ListUtils::create<String>("mzML"));
  }

  ExitCodes main_(int, const char **) override
  {
    int mass_count = getIntOption_("in");
    double min_mass = getDoubleOption_("min_mass");
    double max_mass = getDoubleOption_("max_mass");
    double min_mz = getDoubleOption_("min_mz");
    double max_mz = getDoubleOption_("max_mz");
    double snr = getDoubleOption_("snr");
    String out = getStringOption_("out");



    return EXECUTION_OK;
  }
};

int main(int argc, const char **argv)
{
  TOPPFLASHDeconvSpectrumSimulator tool;
  return tool.main(argc, argv);
}

/// @endcond
