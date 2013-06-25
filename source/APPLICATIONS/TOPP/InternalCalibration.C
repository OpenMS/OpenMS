// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_InternalCalibration InternalCalibration

  @brief Performs an internal calibration on an MS experiment.
 
  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ InternalCalibration \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on MS peak data @n (in mzML format) or feature data </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
  </tr>
  </table>
  </CENTER>
 
  This a simple calibration method: given a list of reference masses and an MS experiment or a feature map,
  the relative errors of the peaks in the data are approximated by linear regression and
  subtracted from the data. The user can choose whether the calibration function shall be
  calculated for each spectrum separately or once for the whole map.
  If this is done scanwise, at least two reference masses need to
  be present in each scan to calculate the calibration function,
  otherwise the spectrum can't be calibrated.
  For the global calibration it is also possible to use a list of (significant) peptide identifications.

  @note The tool assumes the input data is already picked or feature maps.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_InternalCalibration.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_InternalCalibration.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPInternalCalibration :
  public TOPPBase
{
public:
  TOPPInternalCalibration() :
    TOPPBase("InternalCalibration", "Applies an internal calibration.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input peak file ");
    setValidFormats_("in", StringList::create("mzML,featureXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", StringList::create("mzML,featureXML"));
    registerInputFile_("ref_peaks", "<file>", "", "input file containing reference m/z values (either as textfile with one m/z per line and no header or as idXML file)", false);
    setValidFormats_("ref_peaks", StringList::create("csv,idXML"));
    registerStringOption_("type", "<calibration type>", "spectrumwise", "The kind of internal calibration that should be applied.", false);
    setValidStrings_("type", StringList::create("spectrumwise,global"));
    registerOutputFile_("trafo", "<file>", "", "output transformation file (only for global calibration)", false);
    setValidFormats_("trafo", StringList::create("trafoXML"));
    addEmptyLine_();
    registerSubsection_("algorithm", "Settings for the internal calibration.");
  }

  Param getSubsectionDefaults_(const String & /* section*/) const
  {
    Param tmp;
    tmp.insert("", InternalCalibration().getDefaults());
    return tmp;
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String ref = getStringOption_("ref_peaks");
    String type = getStringOption_("type");
    String trafo = getStringOption_("trafo");
    //-------------------------------------------------------------
    // init InternalCalibration
    //-------------------------------------------------------------

    InternalCalibration calib;
    Param param = getParam_().copy("algorithm:", true);
    calib.setParameters(param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    // get reference m/z values
    std::vector<PeptideIdentification> pep_ids;
    vector<DoubleReal> ref_masses;
    bool ids(false);
    if (ref != "")
    {
      ids = FileHandler().getTypeByContent(ref) == FileTypes::IDXML;
      if (ids)
      {
        std::vector<ProteinIdentification> prot_ids;
        IdXMLFile().load(ref, prot_ids, pep_ids);
      }
      else
      {
        TextFile ref_file;
        ref_file.load(ref, true);
        for (TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
        {
          ref_masses.push_back(String(iter->c_str()).toDouble());
        }
      }
    }

    bool features = FileHandler().getTypeByContent(in) == FileTypes::FEATUREXML;
    if (ref == "" && !features)
    {
      std::cout << "Need a file containing the reference peaks!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (type == "spectrumwise" && features)
    {
      std::cout << "Can't perform a spectrumwise calibration on a feature map!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    if (features)
    {
      FeatureMap<> feature_map, calibrated_feature_map;
      FeatureXMLFile f_file;
      f_file.load(in, feature_map);
      if (ref == "")
      {
        std::cout << "Using the peptide identifications stored in the feature map as reference peaks.\n";
        calib.calibrateMapGlobally(feature_map, calibrated_feature_map, trafo);
      }
      else
      {
        std::cout << "Using peptide identifications given with -ref_peaks as reference peaks.\n";
        calib.calibrateMapGlobally(feature_map, calibrated_feature_map, pep_ids, trafo);
      }
      addDataProcessing_(calibrated_feature_map, getProcessingInfo_(DataProcessing::CALIBRATION));
      f_file.store(out, calibrated_feature_map);
      return EXECUTION_OK;
    }


    MSExperiment<Peak1D> ms_exp_raw, ms_exp_calibrated;
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.load(in, ms_exp_raw);

    //-------------------------------------------------------------
    // perform calibration
    //-------------------------------------------------------------
    if (type == "spectrumwise") calib.calibrateMapSpectrumwise(ms_exp_raw, ms_exp_calibrated, ref_masses);
    else if (ids) calib.calibrateMapGlobally(ms_exp_raw, ms_exp_calibrated, pep_ids, trafo);
    else calib.calibrateMapGlobally(ms_exp_raw, ms_exp_calibrated, ref_masses, trafo);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(ms_exp_calibrated, getProcessingInfo_(DataProcessing::CALIBRATION));

    mz_data_file.store(out, ms_exp_calibrated);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPInternalCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
