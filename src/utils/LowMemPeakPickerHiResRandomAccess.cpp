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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

using namespace OpenMS;
using namespace std;

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_LowMemPeakPickerHiResRandomAccess LowMemPeakPickerHiResRandomAccess

  @brief A tool for peak detection in profile data. Executes the peak picking with @ref OpenMS::PeakPickerHiRes "high_res" algorithm.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ PeakPickerHiRes \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> any tool operating on MS peak data @n (in mzML format)</td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterGaussian </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay </td>
  </tr>
  </table>
  </center>

  Reference:\n
  Weisser <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

  The conversion of the "raw" ion count data acquired
  by the machine into peak lists for further processing
  is usually called peak picking or centroiding. The choice of the algorithm
  should mainly depend on the resolution of the data.
  As the name implies, the @ref OpenMS::PeakPickerHiRes "high_res"
  algorithm is fit for high resolution (orbitrap or FTICR) data.

  @ref TOPP_example_signalprocessing_parameters is explained in the TOPP tutorial.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PeakPickerHiRes.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PeakPickerHiRes.html

  For the parameters of the algorithm section see the algorithm documentation: @ref OpenMS::PeakPickerHiRes "PeakPickerHiRes"

  In the following table you, can find example values of the most important algorithm parameters for
  different instrument types. @n These parameters are not valid for all instruments of that type,
  but can be used as a starting point for finding suitable parameters.
  <table>
  <tr BGCOLOR="#EBEBEB">
  <td>&nbsp;</td>
  <td><b>Q-TOF</b></td>
  <td><b>LTQ Orbitrap</b></td>
  </tr>
  <tr>
  <td BGCOLOR="#EBEBEB"><b>signal_to_noise</b></td>
  <td>2</td>
  <td>0</td>
  </tr>
  </table>
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
//
class TOPPLowMemPeakPickerHiResRandomAccess :
  public TOPPBase
{
public:
  TOPPLowMemPeakPickerHiResRandomAccess() :
    TOPPBase("LowMemPeakPickerHiResRandomAccess", "Finds mass spectrometric peaks in profile mass spectra.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input profile data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output peak file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return PeakPickerHiRes().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    IndexedMzMLFileLoader loader; 

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    OnDiscPeakMap exp;
    loader.load(in, exp);

    // We could write out this warning in the constructor if no spectra have come our way ...
    /*
    if (ms_exp_raw.empty() && ms_exp_raw.getChromatograms().size() == 0)
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    */

    // We could check with the first spectrum that we process whether it fulfills the requirements
    // check for peak type (profile data required)
    /*
    if (!ms_exp_raw.empty() && PeakTypeEstimator().estimateType(ms_exp_raw[0].begin(), ms_exp_raw[0].end()) == SpectrumSettings::CENTROID)
    {
      writeLog_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
    }
    */

    // We can check each spectrum and throw an error if it is not sorted
    // check if spectra are sorted
    /*
    for (Size i = 0; i < ms_exp_raw.size(); ++i)
    {
      if (!ms_exp_raw[i].isSorted())
      {
        writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
    */

    // We can check each chromatogram and throw an error if it is not sorted
    //check if chromatograms are sorted
    /*
    for (Size i = 0; i < ms_exp_raw.getChromatograms().size(); ++i)
    {
      if (!ms_exp_raw.getChromatogram(i).isSorted())
      {
        writeLog_("Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
    */

    //-------------------------------------------------------------
    // pick
    //-------------------------------------------------------------
    PeakMap ms_exp_peaks;

    ///////////////////////////////////
    // Create PeakPickerHiRes and hand it to the PPHiResMzMLConsumer
    ///////////////////////////////////
    Param pepi_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to LowMemPeakPickerHiRes", pepi_param, 3);

    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pepi_param);
    pp.pickExperiment(exp, ms_exp_peaks);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    //annotate output with data processing info
    addDataProcessing_(ms_exp_peaks, getProcessingInfo_(DataProcessing::PEAK_PICKING));

    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.store(out, ms_exp_peaks);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPLowMemPeakPickerHiResRandomAccess tool;
  return tool.main(argc, argv);
}

/// @endcond
