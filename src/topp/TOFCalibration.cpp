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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_TOFCalibration TOFCalibration

  @brief Performs an external calibration for tof spectra.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ TOFCalibration \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_InternalCalibration </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
  </tr>
  </table>
  </CENTER>

  Given one or more calibrant spectra containing flight times, the instrument's calibration constants and the
  expected masses the quadratic function \f$y_i = a + bx_i + cx_i^2\f$ is fitted, where \f$x_i\f$ is the ith flight time.
  If there are more than one calibrant spectra the coefficients \f$a\f$, \f$b\f$ and \f$c\f$ are averaged. The fitted function is
  then used to convert the flight times of the given experiment to m/z-values.

  You can choose to calibrate picked or raw data. If you use picked data, set the flag peak_data. If you have
  raw data an additional peak picking step for the calibrant spectra is needed, the parameters for the
  peak picker can be set in the ini-file.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TOFCalibration.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_TOFCalibration.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTOFCalibration :
  public TOPPBase
{
public:
  TOPPTOFCalibration() :
    TOPPBase("TOFCalibration", "Applies time of flight calibration.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input peak or raw data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    addEmptyLine_();
    registerInputFile_("ext_calibrants", "<file>", "", "input file containing the external calibrant spectra (peak or raw data)\n");
    setValidFormats_("ext_calibrants", ListUtils::create<String>("mzML"));
    registerInputFile_("ref_masses", "<file>", "", "input file containing reference masses of the external calibrant spectra (one per line)", true);
    setValidFormats_("ref_masses", ListUtils::create<String>("txt"));
    registerInputFile_("tof_const", "<file>", "", "File containing TOF conversion constants."
                                                  " These can be either two or three constants\n"
                                                  "per set, depending on the conversion type. Either one set for all calibrant spectra \n"
                                                  "(tab separated), or one for each spectrum.\n"
                                                  "For a detailed description, please have a look at the doxygen documentation."
                                                  "(one set, tab separated, per line)", true);
    setValidFormats_("tof_const", ListUtils::create<String>("csv"));
    registerFlag_("peak_data", "set this flag, if you have peak data, not raw data (the picking parameters are accessible only from the INI file).");

    registerSubsection_("algorithm", "Algorithm section for peak picking");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the PeakPicker param
    Param tmp;
    tmp.insert("PeakPicker:", PeakPickerCWT().getDefaults());
    return tmp;
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String in_calib = getStringOption_("ext_calibrants");
    String ref = getStringOption_("ref_masses");
    String conv = getStringOption_("tof_const");
    //-------------------------------------------------------------
    // init TOFCalibration
    //-------------------------------------------------------------

    TOFCalibration calib;
    calib.setLogType(log_type_);
    Param param = getParam_().copy("algorithm:", true);
    calib.setParameters(param);
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    PeakMap ms_exp_calib, ms_exp_raw;
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.load(in_calib, ms_exp_calib);
    mz_data_file.load(in, ms_exp_raw);

    vector<double> ref_masses;
    TextFile ref_file;
    ref_file.load(ref, true);

    for (TextFile::ConstIterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
    {
      ref_masses.push_back(String(iter->c_str()).toDouble());
    }
    TextFile const_file;
    const_file.load(conv, true);
    std::vector<String> vec;
    TextFile::ConstIterator iter = const_file.begin();
    iter->split('\t', vec);

    std::vector<double> ml1, ml2, ml3;
    ml1.push_back(String(vec[0].c_str()).toDouble());
    ml2.push_back(String(vec[1].c_str()).toDouble());
    if (vec.size() == 3)
    {
      ml3.push_back(String(vec[2].c_str()).toDouble());
    }
    ++iter;

    for (; iter != const_file.end(); ++iter)
    {
      iter->split('\t', vec);
      ml1.push_back(String(vec[0].c_str()).toDouble());
      ml2.push_back(String(vec[1].c_str()).toDouble());
      if (vec.size() == 3)
      {
        ml3.push_back(String(vec[2].c_str()).toDouble());
      }
    }

    if (ml1.size() != 1 &&  ml1.size() != ms_exp_calib.size())
    {
      writeLog_("Incorrect number of calibration constants given. Aborting!");
      return INPUT_FILE_CORRUPT;
    }
    calib.setML1s(ml1);
    calib.setML2s(ml2);
    if (!ml3.empty()) calib.setML3s(ml3);

    //-------------------------------------------------------------
    // perform calibration
    //-------------------------------------------------------------
    if (getFlag_("peak_data"))
    {
      calib.calibrate(ms_exp_calib, ms_exp_raw, ref_masses);
    }
    else
    {
      calib.pickAndCalibrate(ms_exp_calib, ms_exp_raw, ref_masses);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(ms_exp_raw, getProcessingInfo_(DataProcessing::CALIBRATION));

    mz_data_file.store(out, ms_exp_raw);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPTOFCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
