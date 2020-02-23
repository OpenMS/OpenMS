// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>
#include <QDir>
#include <algorithm>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_MetFragAdapter MetFragAdapter

  @brief Preprocess Feature or ConsensusMaps for usage in MetFrag

    <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ AssayGeneratorMetabo \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref OpenSWATH pipeline </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref Utils_AccurateMassSearch </td>
          </tr>
      </table>
  </CENTER>

  Preprocess Feature or ConsensusMaps for usage in MetFrag

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES

class TOPPMetFragAdapter :
  public TOPPBase
{
public:
  TOPPMetFragAdapter() :
    TOPPBase("MetFragAdapter", "Preprocess f/cXML for MetFrag", false)
    {}

protected:

  String tmp_dir_ = "";
  String java_executable_ = "";
  String metfrag_executable_ = "";
  int java_memory_ = 3500;

  map<String, int> adduct_str_to_ionmode_ =
      {
        {"[M+H]1+" , 1 },
        {"[M+NH4]1+" , 18 },
        {"[M+Na]1+" , 23 },
        {"[M+K]1+" , 39 },
        {"[M+CH3OH+H]1+" , 33 },
        {"[M+ACN+H]1+" , 42 },
        {"[M+ACN+Na]1+" , 64 },
        {"[M+2ACN+H]1+" , 83 },
        {"[M-H]1-" , -1 },
        {"[M+Cl]1-" , 35 },
        {"[M+HCOO]1-" , 45 },
        {"[M+CH3COO]1-" , 59 },
        {"", 0}
      };

  bool runMetFrag_(const MetaboTargetedAssay& mta, double prec_tol, const String& db = "PubChem")
  {
    double prec_mz = mta.potential_rmts[0].getPrecursorMZ();

    stringstream str;
    for (const auto& rmt : mta.potential_rmts)
    {
      str << String(rmt.getProductMZ()) << "_" << String(rmt.getLibraryIntensity()) << ";";
    }
    String peaklist = str.str();

    QStringList process_params; // the actual process is Java, not MS-GF+!
    QString java_mem = "-Xmx" + QString::number(java_memory_) + "m";

    process_params << java_mem
                   << "-jar" << metfrag_executable_.toQString()
                   << "MetFragPeakListReader=de.ipbhalle.metfraglib.peaklistreader.FilteredStringTandemMassPeakListReader"
                   << "PeakListString=" + peaklist.toQString()
                   << "MetFragDatabaseType=" + db.toQString()
                   << "IonizedPrecursorMass=" + String(prec_mz).toQString()
                   << "DatabaseSearchRelativeMassDeviation=" + QString::number(prec_tol)
                   << "FragmentPeakMatchAbsoluteMassDeviation=0.001"
                   << "FragmentPeakMatchRelativeMassDeviation=5"
                   << "PrecursorIonMode=" + QString::number(adduct_str_to_ionmode_[mta.compound_adduct])
                   // TODO where to get charge from? adduct? parse pos/neg mode somewhere?
                   << "IsPositiveIonMode=TRUE"
                   // TODO allow more scores
                   << "MetFragScoreTypes=FragmenterScore"
                   << "MetFragScoreWeights=1.0"
                   // TODO check if can be streamed
                   << "MetFragCandidateWriter=CSV"
                   // TODO check if unique
                   << "SampleName=" + mta.compound_name.toQString()
                   << "ResultsPath=" + tmp_dir_.toQString()
                   // TODO parameterize both
                   << "MaximumTreeDepth=1"
                   << "MetFragPreProcessingCandidateFilter=UnconnectedCompoundFilter";
    //TODO check more params like NumberThreads

    TOPPBase::ExitCodes exit_code = runExternalProcess_(java_executable_.toQString(), process_params);
    if (exit_code != EXECUTION_OK)
    {
      return false;
    }
    return true;
  }

  void preprocessing_(const String& featureinfo,
                       const MSExperiment& spectra,
                       std::vector<FeatureMap>& v_fp,
                       KDTreeFeatureMaps& fp_map_kd,
                       FeatureMapping::FeatureToMs2Indices& feature_mapping)
  {
    // if fileparameter is given and should be not empty
    if (!featureinfo.empty())
    {
      if (File::exists(featureinfo) && !File::empty(featureinfo))
      {
        // read featureXML
        FeatureXMLFile fxml;
        FeatureMap feature_map;
        fxml.load(featureinfo, feature_map);

        unsigned int num_masstrace_filter = getIntOption_("min_num_masstraces");
        double precursor_mz_tol = getDoubleOption_("precursor_recalibration_window");
        double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
        bool ppm_prec = getStringOption_("precursor_recalibration_window_unit") == "ppm";


        // filter feature by number of masstraces
        auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                                [&num_masstrace_filter](const Feature &feat) -> bool
                                {
                                  unsigned int n_masstraces = feat.getMetaValue("num_of_masstraces");
                                  return n_masstraces < num_masstrace_filter;
                                });
        feature_map.erase(map_it, feature_map.end());

        v_fp.push_back(feature_map);
        fp_map_kd.addMaps(v_fp);

        // mapping of MS2 spectra to features
        feature_mapping = FeatureMapping::assignMS2IndexToFeature(spectra,
                                                                  fp_map_kd,
                                                                  precursor_mz_tol,
                                                                  precursor_rt_tol,
                                                                  ppm_prec);
      }
      else
      {
        throw OpenMS::Exception::FileEmpty(__FILE__,
                                           __LINE__,
                                           __FUNCTION__,
                                           "Error: FeatureXML was empty, please provide a valid file.");
      }
    }
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("executable", "<executable>", "", "Metfrag jar", false, false, ListUtils::create<String>("skipexists"));

    registerInputFileList_("in", "<file(s)>", StringList(), "MzML input file(s) used for assay library generation");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_id", "<file(s)>", StringList(), "FeatureXML input file(s) containing identification information (e.g. AccurateMassSearch)");
    setValidFormats_("in_id", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<file>", "", "Preprocessed spectra to evaluate");
    setValidFormats_("out", ListUtils::create<String>("msp"));

    registerFlag_("use_exact_mass", "Use exact mass for fragment annotation", false);

    // preprocessing
    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);
    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.1, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds]", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);
    registerStringOption_("method", "<choice>", "highest_intensity", "Spectrum with the highest precursor intensity or a consensus spectrum ist used for assay library construction (if no fragment annotation is used).",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));
    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);
    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.90, "Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation", false);
    registerDoubleOption_("transition_threshold", "<num>", 0, "Further transitions need at least x% of the maximum intensity (default 10%)", false);
    registerIntOption_("min_num_masstraces", "<num>", 1, "min num mt", false);


    registerTOPPSubsection_("deisotoping", "deisotoping");
    registerFlag_("deisotoping:use_deisotoper", "Use Deisotoper (if no fragment annotation is used)", false);
    registerDoubleOption_("deisotoping:fragment_tolerance", "<num>", 1, "Tolerance used to match isotopic peaks", false);
    registerStringOption_("deisotoping:fragment_unit", "<choice>", "ppm", "Unit of the fragment tolerance", false);
    setValidStrings_("deisotoping:fragment_unit", ListUtils::create<String>("ppm,Da"));
    registerIntOption_("deisotoping:min_charge", "<num>", 1, "The minimum charge considered", false);
    setMinInt_("deisotoping:min_charge", 1);
    registerIntOption_("deisotoping:max_charge", "<num>", 1, "The maximum charge considered", false);
    setMinInt_("deisotoping:max_charge", 1);
    registerIntOption_("deisotoping:min_isopeaks", "<num>", 2, "The minimum number of isotopic peaks (at least 2) required for an isotopic cluster", false);
    setMinInt_("deisotoping:min_isopeaks", 2);
    registerIntOption_("deisotoping:max_isopeaks", "<num>", 3, "The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster", false);
    setMinInt_("deisotoping:max_isopeaks", 3);
    registerFlag_("deisotoping:keep_only_deisotoped", "Only monoisotopic peaks of fragments with isotopic pattern are retained", false);
    registerFlag_("deisotoping:annotate_charge", "Annotate the charge to the peaks", false);

    registerInputFile_("java_executable", "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, {"is_executable"});
    registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    registerStringOption_("out_workspace_directory", "<directory>", "", "Output directory for MetFrag workspace", false);

  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------
    tmp_dir_ = makeAutoRemoveTempDirectory_();
    java_memory_ = getIntOption_("java_memory");
    java_executable_ = getStringOption_("java_executable");
    metfrag_executable_ = getStringOption_("executable");

    if (!getFlag_("force"))
    {
      if (!JavaInfo::canRun(java_executable_))
      {
        writeLog_("Fatal error: Java is needed to run MS-GF+!");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
    else
    {
      writeLog_("The installation of Java was not checked.");
    }

    // param AssayGeneratorMetabo
    StringList in = getStringList_("in");
    StringList id = getStringList_("in_id");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    bool method_consensus_spectrum = method == "consensus_spectrum";
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");


    double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
    double pre_recal_win = getDoubleOption_("precursor_recalibration_window");
    String pre_recal_win_unit = getStringOption_("precursor_recalibration_window_unit");
    bool ppm_recal = pre_recal_win_unit == "ppm";

    double precursor_mz_distance = getDoubleOption_("precursor_mz_distance");

    double cosine_sim_threshold = getDoubleOption_("cosine_similarity_threshold");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    // param deisotoper
    bool use_deisotoper = getFlag_("deisotoping:use_deisotoper");
    double fragment_tolerance = getDoubleOption_("deisotoping:fragment_tolerance");
    String fragment_unit = getStringOption_("deisotoping:fragment_unit");
    bool fragment_unit_ppm = fragment_unit == "ppm";
    int min_charge = getIntOption_("deisotoping:min_charge");
    int max_charge = getIntOption_("deisotoping:max_charge");
    unsigned int min_isopeaks = getIntOption_("deisotoping:min_isopeaks");
    unsigned int max_isopeaks = getIntOption_("deisotoping:max_isopeaks");
    bool keep_only_deisotoped = getFlag_("deisotoping:keep_only_deisotoped");
    bool annotate_charge = getFlag_("deisotoping:annotate_charge");


    //-------------------------------------------------------------
    // input and check
    //-------------------------------------------------------------

    // check size of .mzML & .featureXML input
    if (in.size() != id.size())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Number of .mzML do not match to the number of .featureXML files. \n Please check and provide the corresponding files.");
    }

    vector<MetaboTargetedAssay> v_mta;

    // iterate over all the files
    for (Size file_counter = 0; file_counter < in.size(); file_counter++)
    {
      // load mzML
      MzMLFile mzml;
      PeakMap spectra;
      mzml.load(in[file_counter], spectra);

      // load featurexml
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(id[file_counter], feature_map);

      // check if featureXML corresponds to mzML
      StringList featurexml_primary_path;
      feature_map.getPrimaryMSRunPath(featurexml_primary_path);

      if (in != featurexml_primary_path)
      {
        OPENMS_LOG_WARN << "Warning: Original paths of the mzML files do not correspond to the featureXML files. Please check and provide the corresponding files." << std::endl;

        OPENMS_LOG_WARN << "Input MzML: " << std::endl;
                        for (const String& it_mzml : in)
                            {
                              OPENMS_LOG_WARN << " " << it_mzml << std::endl;
                            }
        OPENMS_LOG_WARN << "Input FeatureXML: " << std::endl;
                        for (const String& it_fxml : id)
                            {
                              OPENMS_LOG_WARN << " " << it_fxml << std::endl;
                            }
        OPENMS_LOG_WARN << "Original paths: " << std::endl;
                        for (const String& it_fpp : featurexml_primary_path)
                            {
                              OPENMS_LOG_WARN << " " << it_fpp << std::endl;
                            }
      }

      // determine type of spectral data (profile or centroided)
      if (!spectra[0].empty())
      {
        SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

        if (spectrum_type == SpectrumSettings::PROFILE)
        {
          if (!getFlag_("force"))
          {
            throw OpenMS::Exception::FileEmpty(__FILE__,
                                               __LINE__,
                                               __FUNCTION__,
                                               "Error: Profile data provided but centroided spectra expected. ");
          }
        }
      }

      //-------------------------------------------------------------
      // Processing
      //-------------------------------------------------------------

      // sort spectra
      spectra.sortSpectra();

      // check if correct featureXML is given and set use_known_unkowns parameter if no id information is available
      const std::vector<DataProcessing> &processing = feature_map.getDataProcessing();
      for (auto it = processing.begin(); it != processing.end(); ++it)
      {
        if (it->getSoftware().getName() == "FeatureFinderMetabo")
        {
          // if id information is missing set use_known_unknowns to true
          if (feature_map.getProteinIdentifications().empty())
          {
            use_known_unknowns = true;
            OPENMS_LOG_INFO << "Due to the use of data without previous identification "
                     << "use_known_unknowns will be switched on." << std::endl;
          }
        }
      }

      // annotate and recalibrate precursor mz and intensity
      vector<double> delta_mzs;
      vector<double> mzs;
      vector<double> rts;
      PrecursorCorrection::correctToHighestIntensityMS1Peak(spectra, pre_recal_win, ppm_recal, delta_mzs, mzs, rts);

      // always use preprocessing: 
      // run masstrace filter and feature mapping
      vector<FeatureMap> v_fp; // copy FeatureMap via push_back
      KDTreeFeatureMaps fp_map_kd; // reference to *basefeature in vector<FeatureMap>
      FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in vector<FeatureMap>
      preprocessing_(id[file_counter],
                     spectra,
                     v_fp,
                     fp_map_kd,
                     feature_mapping);
    
      // filter known_unknowns based on description (UNKNOWN) (AMS)
      Map<const BaseFeature*, std::vector<size_t>> feature_ms2_spectra_map = feature_mapping.assignedMS2;
      Map<const BaseFeature*, std::vector<size_t>> known_features;
      if (!use_known_unknowns)
      {
        for (auto it = feature_ms2_spectra_map.begin(); it != feature_ms2_spectra_map.end(); ++it)
        {
          const BaseFeature *feature = it->first;
          if (!(feature->getPeptideIdentifications().empty()) &&
              !(feature->getPeptideIdentifications()[0].getHits().empty()))
              {
                String description;
                // one hit is enough for prefiltering
                description = feature->getPeptideIdentifications()[0].getHits()[0].getMetaValue("description");
                // change format of description [name] to name
                description.erase(remove_if(begin(description),
                                            end(description),
                                            [](char c) { return c == '[' || c == ']'; }), end(description));
                known_features.insert({it->first, it->second});
              }
        }
        feature_mapping.assignedMS2 = known_features;
      }

      vector< MetaboTargetedAssay::CompoundSpectrumPair > v_cmp_spec;

      if (use_deisotoper)
      {
        bool make_single_charged = false;
        for (auto& peakmap_it : spectra)
        {
          MSSpectrum& spectrum = peakmap_it;
          if (spectrum.getMSLevel() == 1)
          {
            continue;
          }
          else
          {
            Deisotoper::deisotopeAndSingleCharge(spectrum,
                                                fragment_tolerance,
                                                fragment_unit_ppm,
                                                min_charge,
                                                max_charge,
                                                keep_only_deisotoped,
                                                min_isopeaks,
                                                max_isopeaks,
                                                make_single_charged,
                                                annotate_charge);
          }
        }
      }

      // remove peaks from MS2 which are at a higher mz than the precursor + 10 ppm
      for (auto& peakmap_it : spectra)
      {
        MSSpectrum& spectrum = peakmap_it;
        if (spectrum.getMSLevel() == 1)
        {
          continue;
        }
        else
        {
          // if peak mz higher than precursor mz set intensity to zero
          double prec_mz = spectrum.getPrecursors()[0].getMZ();
          double mass_diff = Math::ppmToMass(10.0, prec_mz);
          for (auto& spec : spectrum)
          {
            if (spec.getMZ() > prec_mz + mass_diff)
            {
              spec.setIntensity(0);
            }
          }
          spectrum.erase(remove_if(spectrum.begin(),
                                   spectrum.end(),
                                   InIntensityRange<PeakMap::PeakType>(1,
                                                                       numeric_limits<PeakMap::PeakType::IntensityType>::max(),
                                                                       true)), spectrum.end());
        }
      }

      // potential transitions of one file
      vector<MetaboTargetedAssay> tmp_mta;
      tmp_mta = MetaboTargetedAssay::extractMetaboTargetedAssay(spectra,
                                                                feature_mapping,
                                                                precursor_rt_tol,
                                                                precursor_mz_distance,
                                                                cosine_sim_threshold,
                                                                transition_threshold,
                                                                method_consensus_spectrum,
                                                                exclude_ms2_precursor,
                                                                file_counter);
      
      // append potential transitions of one file to vector of all files
      v_mta.insert(v_mta.end(), tmp_mta.begin(), tmp_mta.end());
      
    } // end iteration over all files
    //TODO we probably only need one file

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    //MSPFile m;
    //m.store(out, v_mta);

    for (const auto& mta : v_mta)
    {
      //runMetFrag_(mta, precursor_mz_distance);
      runMetFrag_(mta, 5);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPMetFragAdapter tool;
  return tool.main(argc, argv);
}
/// @endcond
