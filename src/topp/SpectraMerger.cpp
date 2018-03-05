// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Chris Bielow, Andreas Bertsch, Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SpectraMerger SpectraMerger

  @brief Allows to add up several spectra.
 
  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SpectraMerger \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
  </tr>
  </table>
  </center>
 
  @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

  This tool can add several consecutive scans, increasing S/N ratio (for MS1 and above) or merge scans which stem from similar precursors (for MS2 and above).

  In any case, the number of scans will be reduced.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SpectraMerger.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SpectraMerger.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraMerger :
  public TOPPBase
{
public:
  TOPPSpectraMerger() :
    TOPPBase("SpectraMerger", "Merges spectra (each MS level separately), increasing S/N ratios.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output mzML file with merged spectra.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerStringOption_("merging_method", "<method>", "average_gaussian", "Method of merging which should be used.", false);
    setValidStrings_("merging_method", ListUtils::create<String>("average_gaussian,average_tophat,precursor_method,block_method"));

    registerSubsection_("algorithm", "Algorithm section for merging spectra");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return SpectraMerger().getParameters();
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String out(getStringOption_("out"));
    String merging_method(getStringOption_("merging_method"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    FileHandler fh;
    FileTypes::Type in_type = fh.getType(in);

    PeakMap exp;
    fh.loadExperiment(in, exp, in_type, log_type_);
    exp.sortSpectra();

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    SpectraMerger merger;
    merger.setLogType(log_type_);
    merger.setParameters(getParam_().copy("algorithm:", true));
    if (merging_method == "precursor_method")
    {
      merger.mergeSpectraPrecursors(exp);
    }
    else if (merging_method == "block_method")
    {
      merger.mergeSpectraBlockWise(exp);
    }
    else if (merging_method == "average_gaussian")
    {
      merger.average(exp, "gaussian");
    }
    else if (merging_method == "average_tophat")
    {
      merger.average(exp, "tophat");
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    fh.storeExperiment(out, exp, log_type_);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPSpectraMerger tool;
  return tool.main(argc, argv);
}

/// @endcond
