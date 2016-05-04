// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

  Given reference masses (as either peptide identifications or as list of fixed masses) an MS experiment or a feature map
  can be recalibrated using a linear regression.
  
  For MS experiments, the user can choose whether the calibration function shall be
  calculated for each spectrum separately or once for the whole map.
  If this is done scan-wise, at least two reference masses need to
  be present in each scan to calculate the calibration function.

  Feature maps must be calibrated globally.
  They can only be calibrated using either external idXML or their internal peptide annotation.
  Make sure that each feature has a single peptide identification attached to it (see @ref TOPP_IDFilter) to make the 
  reference mass unambiguous.


  @note The tool assumes the input data is already picked or feature maps.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

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
    // data
    registerInputFile_("in", "<file>", "", "Input peak file");
    setValidFormats_("in", ListUtils::create<String>("mzML,featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML,featureXML"));
    
    // transformation
    registerInputFile_("ref_peaks", "<file>", "", "Input file containing reference m/z values (either as text file with one m/z per line and no header or as idXML file)", false);
    setValidFormats_("ref_peaks", ListUtils::create<String>("csv,idXML"));
    
       
    registerStringOption_("scope", "<calibration type>", "spectrum", "The kind of internal calibration that should be applied.", false);
    setValidStrings_("scope", ListUtils::create<String>("spectrum,global"));

    registerOutputFile_("out_trafo", "<file>", "", "Output transformation file (only for global calibration)", false);
    setValidFormats_("out_trafo", ListUtils::create<String>("trafoXML"));
    addEmptyLine_();
    registerSubsection_("algorithm", "Settings for the internal calibration.");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    Param tmp;
    tmp.insert("", InternalCalibration().getDefaults());
    return tmp;
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String ref = getStringOption_("ref_peaks");
    String scope = getStringOption_("scope");
    String out_trafo = getStringOption_("out_trafo");
    
    //-------------------------------------------------------------
    // init InternalCalibration
    //-------------------------------------------------------------
    InternalCalibration calib;
    Param param = getParam_().copy("algorithm:", true);
    calib.setParameters(param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    // featureXML input
    if (FileHandler().getTypeByContent(in) == FileTypes::FEATUREXML)
    {
      if (scope == "spectrum")
      {
        LOG_ERROR << "Can't perform a spectrum-wise calibration on a feature map!" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      FeatureMap feature_map;
      FeatureXMLFile f_file;
      f_file.load(in, feature_map);
      if (ref == "")
      {
        std::cout << "Using the peptide identifications stored in the feature map as reference peaks.\n";
        calib.calibrateMapGlobally(feature_map, out_trafo);
      }
      else
      {
        std::cout << "Using peptide identifications given with -ref_peaks as reference peaks.\n";
        std::vector<PeptideIdentification> pep_ids;
        if (FileHandler().getTypeByContent(ref) == FileTypes::IDXML)
        {
          std::vector<ProteinIdentification> prot_ids;
          IdXMLFile().load(ref, prot_ids, pep_ids);
        } else
        {
          LOG_ERROR << "Input file -ref_peaks is not of type idXML. Please provide idXML as input!" << std::endl;
          return ILLEGAL_PARAMETERS;
        }
        calib.calibrateMapGlobally(feature_map, pep_ids, out_trafo);
      }
      addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::CALIBRATION));
      f_file.store(out, feature_map);
      return EXECUTION_OK;
    }

    // Raw data
    MSExperiment<Peak1D> ms_exp_raw;
    MzMLFile mz_file;
    mz_file.setLogType(log_type_);
    mz_file.load(in, ms_exp_raw);

    if (ref.empty())
    {
      std::cout << "Need a file containing the reference peaks provided via -ref_peaks input argument!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    if (FileHandler().getTypeByContent(ref) == FileTypes::IDXML)
    {
      std::vector<PeptideIdentification> pep_ids;
      std::vector<ProteinIdentification> prot_ids;
      IdXMLFile().load(ref, prot_ids, pep_ids);
      calib.calibrateMapGlobally(ms_exp_raw, pep_ids, out_trafo);
    }
    else
    {
      TextFile ref_file;
      ref_file.load(ref, true);
      vector<double> ref_masses;
      for (TextFile::ConstIterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
      {
        ref_masses.push_back(String(iter->c_str()).toDouble());
      }
      if (scope == "spectrum") calib.calibrateMapSpectrumwise(ms_exp_raw, ref_masses);
      else calib.calibrateMapGlobally(ms_exp_raw, ref_masses, out_trafo);
    }


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(ms_exp_raw, getProcessingInfo_(DataProcessing::CALIBRATION));

    mz_file.store(out, ms_exp_raw);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPInternalCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
