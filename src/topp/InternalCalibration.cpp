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
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>

#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_InternalCalibration InternalCalibration

  @brief Performs an internal mass recalibration on an MS experiment.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ InternalCalibration \f$ \longrightarrow \f$</td>
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

  Given reference masses (as either peptide identifications or as list of fixed masses) an MS experiment
  can be recalibrated using a linear or quadratic regression fitted to the observed vs. the theoretical masses.

  Chose one of two optional input files:
   1) peptide identifications (from featureXML or idXML) using 'id_in'
   2) lock masses using 'lock_in'
  
  The user can choose whether the calibration function shall be
  calculated for each spectrum separately or once for the whole map.
  If this is done scan-wise, a user-defined range of neighboring spectra
  is searched for lock masses/peptide IDs. They are used to build a model, which is applied to
  the spectrum at hand.
  The RT range ('RT_chunking') should be small enough to resolve time-dependent change of decalibration, but wide enough
  to have enough calibrant masses for a stable model. A linear model requires at least two calibrants, a quadradic at least three.
  Usually, the RT range should provide about 3x more calibrants than required, i.e. 6(=3x2) for linear, and 9(=3x3) for quadratic models.
  If the calibrant data is too sparse for a certain scan, the closest neighboring model will be used automatically.
  If no model can be calculated anywhere, the tool will fail.
  
  Optional quality control output files allow to judge the success of calibration. It is strongly advised to inspect them.
  If PNG images are requested, 'R' (statistical programming language) needs to be installed and available on the system path!

  Outlier detection is supported using the RANSAC algorithm. However, usually it's better to provide high-confidence calibrants instead of
  relying on automatic removal of outliers.

  Post calibration statistics (median ppm and median-absolute-deviation) are automatically computed.
  The calibration is deemed successful if the statistics are within certain bounds ('goodness:XXX').


  Detailed description for each calibration method:
  1) [id_in] The peptide identifications should be derived from the very same mzML file using a wide precursor window (e.g. 25 ppm), which captures
     the possible decalibration. Subsequently, the IDs should be filtered for high confidence (e.g. low FDR, ideally FDR=0.0) and given as input to this tool.
     Remaining outliers can be removed by using RANSAC.
     The data might benefit from a precursor mass correction (e.g. using @ref TOPP_HighResPrecursorMassCorrector), before an MS/MS search is done.
     The list of calibrants is derived solely from the idXML/featureXML and only the resulting model is applied to the mzML.
  
  2) [lock_in] Calibration can be performed using specific lock masses which occur in most spectra. The structure of the cal:lock_in CSV file is as follows:
    Each line represents one lock mass in the format: <m/z>, <ms-level>, <charge>
    Lines starting with # are treated as comments and ignored. The ms-level is usually '1', but you can also use '2' if there are fragment ions commonly occurring.

    Example:
    @code
      # lock mass at 574 m/z at MS1 with charge 2
      574.345, 1, 2
    @endcode

    Additional filters ('cal:lock_require_mono', 'cal:lock_require_iso') allow to exclude spurious false-positive calibrant peaks.
    These filters require knowledge of the charge state, thus charge needs to be specified in the input CSV.
    Detailed information on which lock masses passed these filters are available when -debug is used (any level).

    The calibration function will use all lock masses (i.e. from all ms-levels) within the defined RT range to calibrate a spectrum. Thus, care should be taken that
    spectra from ms-levels specified here, are recorded using the same mass analyzer (MA). This is no issue for a Q-Exactive (which only has one MA),
    but depends on the acquisition scheme for instruments with two/three MAs (e.g. for Orbitrap Velos, MS/MS spectra are commonly acquired in the ion trap and should not be used during calibration of MS1).

  General remarks:
  The user can select what MS levels are subjected to calibration. Calibration must be done once for each mass analyzer.
  Usually, peptide ID's provide calibration points for MS1 precursors, i.e. are suitable for MS1. They are applicable for MS2 only if
  the same mass analyzer was used (e.g. Q-Exactive). In other words, MS/MS spectra acquired using the ion trap analyzer of a Velos cannot be calibrated using
  peptide ID's.
  Precursor m/z associated to higher-level MS spectra are corrected if their precursor spectra are subject to calibration, 
  e.g. precursor information within MS2 spectra is calibrated if target ms-level is set to 1.
  Lock masses ('cal:lock_in') can be specified freely for MS1 and/or MS2.


  @note The tool assumes the input data is already picked/centroided.

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
    TOPPBase("InternalCalibration", "Applies an internal mass recalibration.")
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
    registerInputFile_("rscript_executable", "<file>", "Rscript", "Path to the Rscript executable (default: 'Rscript').", false);
        
    addEmptyLine_();

    registerDoubleOption_("ppm_match_tolerance", "<delta m/z in [ppm]>", 25, "Finding calibrants in raw data uses this tolerance (for lock masses and ID's).", false);

    // transformation
    registerTOPPSubsection_("cal", "Chose one of two optional input files ('id_in' or 'lock_in') to define the calibration masses/function");
    registerInputFile_("cal:id_in", "<file>", "", "Identifications or features whose peptide ID's serve as calibration masses.", false);
    setValidFormats_("cal:id_in", ListUtils::create<String>("idXML,featureXML"));
    registerInputFile_("cal:lock_in", "<file>", "", "Input file containing reference m/z values (text file with each line as: m/z ms-level charge) which occur in all scans.", false);
    setValidFormats_("cal:lock_in", ListUtils::create<String>("csv"));
    registerOutputFile_("cal:lock_out", "<file>", "", "Optional output file containing peaks from 'in' which were matched to reference m/z values. Useful to see which peaks were used for calibration.", false);
    setValidFormats_("cal:lock_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("cal:lock_fail_out", "<file>", "", "Optional output file containing lock masses which were NOT found or accepted(!) in data from 'in'. Useful to see which peaks were used for calibration.", false);
    setValidFormats_("cal:lock_fail_out", ListUtils::create<String>("mzML"));
    registerFlag_("cal:lock_require_mono", "Require all lock masses to be monoisotopic, i.e. not the iso1, iso2 etc ('charge' column is used to determine the spacing). Peaks which are not mono-isotopic are not used.");
    registerFlag_("cal:lock_require_iso", "Require all lock masses to have at least the +1 isotope. Peaks without isotope pattern are not used.");
    registerStringOption_("cal:model_type", 
                          "<model>",
                          MZTrafoModel::enumToName(MZTrafoModel::LINEAR_WEIGHTED),
                          "Type of function to be fitted to the calibration points.",
                          false);
    setValidStrings_("cal:model_type", MZTrafoModel::names_of_modeltype, MZTrafoModel::SIZE_OF_MODELTYPE);

    addEmptyLine_();
    
    registerIntList_("ms_level", "i j ...", ListUtils::create<int>("1,2,3"), "Target MS levels to apply the transformation onto. Does not affect calibrant collection.", false);
    
    registerDoubleOption_("RT_chunking", "<RT window in [sec]>", 300, "RT window (one-sided, i.e. left->center, or center->right) around an MS scan in which calibrants are collected to build a model. Set to -1 to use ALL calibrants for all scans, i.e. a global model.", false);
    
    registerTOPPSubsection_("RANSAC", "Robust outlier removal using RANSAC");
    registerFlag_("RANSAC:enabled", "Apply RANSAC to calibration points to remove outliers before fitting a model.");
    // RANSAC:n is automatically taken from the input model (i.e. n=2 for linear, n=3 for quadratic)
    //registerIntOption_("RANSAC:pc_n", "<# points>", 20, "Percentage (1-99) of initial model points from available data.", false);
    //setMinInt_("RANSAC:pc_n", 1);
    //setMaxInt_("RANSAC:pc_n", 99);
    registerDoubleOption_("RANSAC:threshold", "<threshold>", 10.0, "Threshold for accepting inliers (instrument precision (not accuracy!) as ppm^2 distance)", false);
    registerIntOption_("RANSAC:pc_inliers", "<# inliers>", 30, "Minimum percentage (of available data) of inliers (<threshold away from model) to accept the model.", false);
    setMinInt_("RANSAC:pc_inliers", 1);
    setMaxInt_("RANSAC:pc_inliers", 99);
    registerIntOption_("RANSAC:iter", "<# iterations>", 70, "Maximal # iterations.", false);
    
    registerTOPPSubsection_("goodness", "Thresholds for accepting calibration success");
    registerDoubleOption_("goodness:median", "<threshold>", 4.0, "The median ppm error of calibrated masses must be smaller than this threshold.", false);
    registerDoubleOption_("goodness:MAD", "<threshold>", 1.0, "The median absolute deviation of the ppm error of calibrated masses must be smaller than this threshold.", false);

    registerTOPPSubsection_("quality_control", "Tables and plots to verify calibration performance");
    registerOutputFile_("quality_control:models", "<table>", "", "Table of model parameters for each spectrum.", false);
    setValidFormats_("quality_control:models", ListUtils::create<String>("csv"));
    registerOutputFile_("quality_control:models_plot", "<image>", "", "Plot image of model parameters for each spectrum.", false);
    setValidFormats_("quality_control:models_plot", ListUtils::create<String>("png"));
    registerOutputFile_("quality_control:residuals", "<table>", "", "Table of pre- and post calibration errors.", false);
    setValidFormats_("quality_control:residuals", ListUtils::create<String>("csv"));
    registerOutputFile_("quality_control:residuals_plot", "<image>", "", "Plot image of pre- and post calibration errors.", false);
    setValidFormats_("quality_control:residuals_plot", ListUtils::create<String>("png"));


  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String rscript_executable = getStringOption_("rscript_executable"); 
    String cal_id = getStringOption_("cal:id_in");
    String cal_lock = getStringOption_("cal:lock_in");
    String file_cal_lock_out = getStringOption_("cal:lock_out");
    String file_cal_lock_fail_out = getStringOption_("cal:lock_fail_out");
    double rt_chunk = getDoubleOption_("RT_chunking");
    
    IntList ms_level = getIntList_("ms_level");

    if (((int)!cal_lock.empty() + (int)!cal_id.empty()) != 1)
    {
      LOG_ERROR << "Conflicting input given. Please provide only ONE of either 'cal:id_in' or 'cal:lock_in'!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    // Raw data
    PeakMap exp;
    MzMLFile mz_file;
    mz_file.setLogType(log_type_);
    mz_file.load(in, exp);

    InternalCalibration ic;
    ic.setLogType(log_type_);

    double tol_ppm = getDoubleOption_("ppm_match_tolerance");

    // featureXML/idXML input
    if (!cal_id.empty())
    {
      FileTypes::Type ftype = FileHandler().getTypeByContent(cal_id);
      if (ftype == FileTypes::FEATUREXML)
      {
        FeatureMap feature_map;
        FeatureXMLFile().load(cal_id, feature_map);
        ic.fillCalibrants(feature_map, tol_ppm);
      }
      else if (ftype == FileTypes::IDXML)
      {
        std::vector<ProteinIdentification> prot_ids;
        std::vector<PeptideIdentification> pep_ids; 
        IdXMLFile().load(cal_id, prot_ids, pep_ids);
        ic.fillCalibrants(pep_ids, tol_ppm);
      }
    }
    else if (!cal_lock.empty())
    { // CSV file calibrant masses
      // load CSV
      TextFile ref_file;
      ref_file.load(cal_lock, true, -1, true);
      vector<InternalCalibration::LockMass> ref_masses;
      for (TextFile::ConstIterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
      {
        if (iter->hasPrefix("#")) continue;
        // each line has:
        //   m/z, ms-level, charge
        std::vector<String> vec;
        iter->split(",", vec);
        if (vec.size() != 3)
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Input file ") + cal_lock + " does not have three comma-separated entries per row!");
        }
        ref_masses.push_back(InternalCalibration::LockMass(vec[0].toDouble(), vec[1].toInt(), vec[2].toInt()));
      }

      bool lock_require_mono = getFlag_("cal:lock_require_mono");
      bool lock_require_iso = getFlag_("cal:lock_require_iso");

      // match calibrants to data
      CalibrationData failed_points;
      ic.fillCalibrants(exp, ref_masses, tol_ppm, lock_require_mono, lock_require_iso, failed_points, debug_level_ > 0);
      
      // write matched lock mass peaks
      if (!file_cal_lock_out.empty())
      {
        LOG_INFO << "\nWriting matched lock masses to mzML file '" << file_cal_lock_out << "'." << std::endl;
        PeakMap exp_out;
        exp_out.set2DData(ic.getCalibrationPoints(), CalibrationData::getMetaValues());
        mz_file.store(file_cal_lock_out, exp_out);
      }
      if (!file_cal_lock_fail_out.empty())
      {
        LOG_INFO << "\nWriting unmatched lock masses to mzML file '" << file_cal_lock_fail_out << "'." << std::endl;
        PeakMap exp_out;
        exp_out.set2DData(failed_points, CalibrationData::getMetaValues());
        mz_file.store(file_cal_lock_fail_out, exp_out);
      }
    }
    
    bool use_RANSAC = getFlag_("RANSAC:enabled");

    if (ic.getCalibrationPoints().empty())
    {
      LOG_ERROR << "No calibration points found! Check your Raw data and calibration masses. Aborting!" << std::endl;
      return UNEXPECTED_RESULT;
    }
      
    //
    // create models and calibrate
    //
    String model_type = getStringOption_("cal:model_type");
    MZTrafoModel::MODELTYPE md = MZTrafoModel::nameToEnum(model_type);
    Size RANSAC_initial_points = model_type.hasSubstring("linear") ? 2 : 3;
    Math::RANSACParam p(RANSAC_initial_points, getIntOption_("RANSAC:iter"), getDoubleOption_("RANSAC:threshold"), getIntOption_("RANSAC:pc_inliers"), true);
    MZTrafoModel::setRANSACParams(p);
    // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
    MZTrafoModel::setCoefficientLimits(tol_ppm, tol_ppm, 0.5); 

    if (!ic.calibrate(exp, ms_level, md, rt_chunk, use_RANSAC, 
                      getDoubleOption_("goodness:median"),
                      getDoubleOption_("goodness:MAD"), 
                      getStringOption_("quality_control:models"),
                      getStringOption_("quality_control:models_plot"),
                      getStringOption_("quality_control:residuals"),
                      getStringOption_("quality_control:residuals_plot"),
                      rscript_executable))
    {
      LOG_ERROR << "\nCalibration failed. See error message above!" << std::endl;
      return UNEXPECTED_RESULT;
    }
 
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
  TOPPInternalCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
