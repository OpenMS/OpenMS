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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_ExternalCalibration ExternalCalibration

  @brief Performs an mass recalibration on an MS experiment using an external calibration function.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ExternalCalibration \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on MS peak data @n (in mzML format)</td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
  </tr>
  </table>
  </CENTER>

  Recalibrates an MS experiment globally using a constant, linear or quadratic calibration on the observed ppm error,
  i.e. using offset=-5, slope=0, power=0 assumes the observed data has -5 ppm decalibration, i.e. the observed m/z is too small and should be increased by 5 ppm! Slope and power
  are coefficients for the observed m/z, i.e. y = offset + x * slope + x * x * power, where x is the observed m/z and y 
  is the resulting mass correction in ppm. Usually slope and offset are very small (< 0.01).
  If you only want a 'rough' recalibration, using offset is usually sufficient.
  
  The calibration function is applied to all scans (with the desired level, see below), i.e. time dependent recalibration cannot be modeled.
    
  The user can select what MS levels are subjected to calibration.
  Spectra with other MS levels remain unchanged.
  Calibration must be done once for each mass analyzer.

  Either raw or centroided data can be used as input.
  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_ExternalCalibration.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_ExternalCalibration.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES



class TOPPExternalCalibration :
  public TOPPBase
{
public:
  TOPPExternalCalibration() :
    TOPPBase("ExternalCalibration", "Applies an external mass recalibration.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    // data
    registerInputFile_("in", "<file>", "", "Input peak file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
        
    addEmptyLine_();

    registerDoubleOption_("offset", "", 0.0, "Mass offset in ppm", false);
    registerDoubleOption_("slope", "", 0.0, "Slope (dependent on m/z)", false);
    registerDoubleOption_("power", "", 0.0, "Power (dependent on m/z)", false);
    
    addEmptyLine_();
    
    registerIntList_("ms_level", "i j ...", ListUtils::create<int>("1,2,3"), "Target MS levels to apply the transformation onto. Scans with other levels remain unchanged.", false);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out"); 
    
    IntList ms_level = getIntList_("ms_level");

    double offset = getDoubleOption_("offset");
    double slope = getDoubleOption_("slope");
    double power = getDoubleOption_("power");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    // Raw data
    PeakMap exp;
    MzMLFile mz_file;
    mz_file.setLogType(log_type_);
    mz_file.load(in, exp);

    MZTrafoModel tm;
    tm.setCoefficients(offset, slope, power);

    InternalCalibration ic;
    ic.setLogType(log_type_);
    ic.applyTransformation(exp, ms_level, tm);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::CALIBRATION));

    mz_file.store(out, exp);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPExternalCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
