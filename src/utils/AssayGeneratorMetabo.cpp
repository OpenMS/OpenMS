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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <map> //insert

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h> 
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_AssayGeneratorMetabo AssayGeneratorMetabo

  @brief Generates an assay library using DDA data (Metabolomics)

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

  Build an assay library from DDA data (MS and MS/MS) (mzML).
  Please provide a list of features found in the data (featureXML).

  Features can be detected using the FeatureFinderMetabo (FFM) and identifcation information
  can be added using the AccurateMassSearch feautreXML output.

  If the FFM featureXML is provided the "use_known_unknowns" flag is used automatically.

  Internal procedure AssayGeneratorMetabo: \n
  1. Input mzML and featureXML \n
  2. Reannotate precursor mz and intensity \n
  3. Filter feature by number of masstraces \n
  4. Assign precursors to specific feature (FeatureMapping) \n
  5. Extract feature meta information (if possible) \n
  6. Find MS2 spectrum with highest intensity precursor for one feature \n
  7. Dependent on the method fragment annotation via SIRIUS is used for transition
  extraction. \n
  If not fragment annotation is performed either the MS2 with the highest intensity precursor or a consensus spectrum
   can be used for the transition extractuib. \n
  8. Calculate thresholds (maximum and minimum intensity for transition peak) \n
  9. Extract and write transitions (tsv, traml) \n

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES

class TOPPAssayGeneratorMetabo :
  public TOPPBase,
  private TransitionTSVFile
{
public:
  TOPPAssayGeneratorMetabo() :
    TOPPBase("AssayGeneratorMetabo", "Assay library generation from DDA data (Metabolomics)", false)
    {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("executable", "<executable>", "", "SIRIUS executable e.g. sirius", false, false, ListUtils::create<String>("skipexists"));

    registerInputFileList_("in", "<file(s)>", StringList(), "MzML input file(s) used for assay library generation");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_id", "<file(s)>", StringList(), "FeatureXML input file(s) containing identification information (e.g. AccurateMassSearch)");
    setValidFormats_("in_id", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv,traML,pqp"));

    registerStringOption_("fragment_annotation", "<choice>", "none", "Fragment annotation method",false);
    setValidStrings_("fragment_annotation", ListUtils::create<String>("none,sirius"));

    registerStringOption_("method", "<choice>", "highest_intensity", "Spectrum with the highest precursor intensity or a consensus spectrum ist used for assay library construction (if no fragment annotation is used).",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerFlag_("use_exact_mass", "Use exact mass for fragment annotation", false);
    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);

    // preprocessing
    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);
    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.1, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds]", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);

    // transition extraction 
    registerIntOption_("min_transitions", "<int>", 3, "minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "maximal number of transitions", false);
    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.98, "Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation", false);
    registerDoubleOption_("transition_threshold", "<num>", 10, "Further transitions need at least x% of the maximum intensity (default 10%)", false);

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

    // sirius 
    registerFullParam_(SiriusAdapterAlgorithm().getDefaults());
  }

  static bool extractAndCompareScanIndexLess_(const String& i, const String& j)
  {
    return (SiriusMzTabWriter::extract_scan_index(i) < SiriusMzTabWriter::extract_scan_index(j));
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    // param AssayGeneratorMetabo
    StringList in = getStringList_("in");
    StringList id = getStringList_("in_id");
    String out = getStringOption_("out");
    String fragment_annotation = getStringOption_("fragment_annotation");
    String method = getStringOption_("method");
    bool use_fragment_annotation = fragment_annotation == "sirius" ? true : false;
    bool method_consensus_spectrum = method == "consensus_spectrum" ? true : false;
    bool use_exact_mass = getFlag_("use_exact_mass");
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");

    int min_transitions = getIntOption_("min_transitions");
    int max_transitions = getIntOption_("max_transitions");

    double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
    double pre_recal_win = getDoubleOption_("precursor_recalibration_window");
    String pre_recal_win_unit = getStringOption_("precursor_recalibration_window_unit");
    bool ppm_recal = pre_recal_win_unit == "ppm" ? true : false;

    double precursor_mz_distance = getDoubleOption_("precursor_mz_distance");

    double cosine_sim_threshold = getDoubleOption_("cosine_similarity_threshold");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    // param deisotoper
    bool use_deisotoper = getFlag_("deisotoping:use_deisotoper");
    double fragment_tolerance = getDoubleOption_("deisotoping:fragment_tolerance");
    String fragment_unit = getStringOption_("deisotoping:fragment_unit");
    bool fragment_unit_ppm = fragment_unit == "ppm" ? true : false;
    int min_charge = getIntOption_("deisotoping:min_charge");
    int max_charge = getIntOption_("deisotoping:max_charge");
    unsigned int min_isopeaks = getIntOption_("deisotoping:min_isopeaks");
    unsigned int max_isopeaks = getIntOption_("deisotoping:max_isopeaks");
    bool keep_only_deisotoped = getFlag_("deisotoping:keep_only_deisotoped");
    bool annotate_charge = getFlag_("deisotoping:annotate_charge");

    // param SiriusAdapterAlgorithm
    String executable = getStringOption_("executable");
    Param combined; 
    SiriusAdapterAlgorithm sirius_algo;
    Param preprocessing = getParam_().copy("preprocessing", false);
    Param sirius = getParam_().copy("sirius", false);
    combined.insert("", preprocessing);
    combined.insert("", sirius);
    sirius_algo.setParameters(combined);

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
    for (unsigned file_counter = 0; file_counter < in.size(); file_counter++)
    {
      // load mzML
      MzMLFile mzml;
      PeakMap spectra;
      mzml.load(in[file_counter], spectra);

      // load featurexml
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(id[file_counter], feature_map);

      // need featureXML with Sourcefile have a look additional to ams
      StringList mzml_primary_path;
      StringList featurexml_primary_path;
      spectra.getPrimaryMSRunPath(mzml_primary_path);
      feature_map.getPrimaryMSRunPath(featurexml_primary_path);

      if (mzml_primary_path != featurexml_primary_path)
      {
        throw Exception::MissingInformation(__FILE__,
                                            __LINE__,
                                            OPENMS_PRETTY_FUNCTION,
                                            "Path of the original input file do not match in the .mzML and .featureXML files. \n Please check and provide the corresponding files.");
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
            LOG_INFO << "Due to the use of data without previous identification "
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
      SiriusAdapterAlgorithm::preprocessingSirius(id[file_counter],
                                                  spectra,
                                                  v_fp,
                                                  fp_map_kd,
                                                  sirius_algo,
                                                  feature_mapping);
    
      // filter known_unkowns based on description (UNKNOWN) (AMS)
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
      if (use_fragment_annotation && executable.empty())
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "SIRIUS executable was not found.");
      }
      else if (use_fragment_annotation && !executable.empty())
      {
        // make temporary files
        SiriusAdapterAlgorithm::SiriusTmpStruct sirius_tmp = SiriusAdapterAlgorithm::constructSiriusTmpStruct();
        String tmp_dir = sirius_tmp.tmp_dir;
        String tmp_ms_file = sirius_tmp.tmp_ms_file;
        String tmp_out_dir = sirius_tmp.tmp_out_dir;
  
        // write msfile and store the compound information in CompoundInfo Object
        vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
        bool feature_only = (sirius_algo.getFeatureOnly() == "true") ? true : false;
        bool no_mt_info = (sirius_algo.getNoMasstraceInfoIsotopePattern() == "true") ? true : false;
        int isotope_pattern_iterations = sirius_algo.getIsotopePatternIterations();
        SiriusMSFile::store(spectra,
                            tmp_ms_file,
                            feature_mapping,
                            feature_only,
                            isotope_pattern_iterations,
                            no_mt_info,
                            v_cmpinfo);

        SiriusAdapterAlgorithm::checkFeatureSpectraNumber(id[file_counter],
                                                          feature_mapping,
                                                          spectra,
                                                          sirius_algo);
  
        // calls SIRIUS and returns vector of paths to sirius folder structure
        std::vector<String> subdirs;
        String out_csifingerid;
        subdirs = SiriusAdapterAlgorithm::callSiriusQProcess(tmp_ms_file,
                                                             tmp_out_dir,
                                                             executable,
                                                             out_csifingerid,
                                                             sirius_algo);
  
        // sort vector path list
        std::sort(subdirs.begin(), subdirs.end(), extractAndCompareScanIndexLess_);
        LOG_DEBUG << subdirs.size() << " spectra were annotated using SIRIUS." << std::endl;
  
        // get Sirius FragmentAnnotion from subdirs
        vector<MSSpectrum> annotated_spectra;
        for (const auto& subdir : subdirs)
        {
          MSSpectrum annotated_spectrum;
          SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(subdir, 
                                                                           annotated_spectrum, 
                                                                           use_exact_mass);
          annotated_spectra.push_back(std::move(annotated_spectrum));
        }
        
        // clean tmp directory if debug level < 2 
        if (debug_level_ >= 2)
        {
          writeDebug_("Keeping temporary files in directory '" + tmp_dir + " and msfile at this location "+ tmp_ms_file + ". Set debug level to 1 or lower to remove them.", 2);
        }
        else
        {
          if (tmp_dir.empty() == false)
          {
            writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 0);
            File::removeDir(tmp_dir.toQString());
          }
          if (tmp_ms_file.empty() == false)
          {
            writeDebug_("Deleting temporary msfile '" + tmp_ms_file + "'. Set debug level to 2 or higher to keep it.", 0);
            File::remove(tmp_ms_file); // remove msfile
          }
        }

        // pair compoundInfo and fragment annotation msspectrum (using the mid)
        for (const auto& cmp : v_cmpinfo)
        {
          for (const auto& spec_fa : annotated_spectra)
          {
            // mid is saved in Name of the spectrum
            if (std::any_of(cmp.mids.begin(), cmp.mids.end(), [spec_fa](const String &str){ return str == spec_fa.getName();}))
            {
              MetaboTargetedAssay::CompoundSpectrumPair csp;
              csp.compoundspectrumpair = std::make_pair(cmp,spec_fa);
              v_cmp_spec.push_back(std::move(csp));
            }
          }
        }
      }
      else // use heuristic
      {
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

        // remove peaks form MS2 which are at a higher mz than the precursor + 10 ppm
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
      }

      // potential transitions of one file
      vector<MetaboTargetedAssay> tmp_mta;
      if (use_fragment_annotation)
      {
        tmp_mta = MetaboTargetedAssay::extractMetaboTargetedAssayFragmentAnnotation(v_cmp_spec,
                                                                                    transition_threshold,
                                                                                    use_exact_mass,
                                                                                    exclude_ms2_precursor,
                                                                                    file_counter);
      }
      else // use heuristics
      {
        tmp_mta = MetaboTargetedAssay::extractMetaboTargetedAssay(spectra,
                                                                  feature_mapping,
                                                                  precursor_rt_tol,
                                                                  precursor_mz_distance,
                                                                  cosine_sim_threshold,
                                                                  transition_threshold,
                                                                  method_consensus_spectrum,
                                                                  exclude_ms2_precursor,
                                                                  file_counter);
      }
      
      // append potential transitions of one file to vector of all files
      v_mta.insert(v_mta.end(), tmp_mta.begin(), tmp_mta.end());
      
    } // end iteration over all files

    // use first rank based on precursor intensity
    std::map< std::pair <String,String>, MetaboTargetedAssay > map_mta;
    for (const auto& it : v_mta)
    {
      pair<String,String> pair_mta = make_pair(it.compound_name, it.compound_adduct);

      // check if value in map with key k does not exists and fill with current pair
      if (map_mta.count(pair_mta) == 0)
      {
        map_mta[pair_mta] = it;
      }
      else
      {
        // check which on has the higher intensity precursor and if replace the current value
        double map_precursor_int = map_mta.at(pair_mta).precursor_int;
        double current_precursor_int = it.precursor_int;
        if (map_precursor_int < current_precursor_int)
        {
          map_mta[pair_mta] = it;
        }
      }
    }

    // merge possible transitions
    vector<TargetedExperiment::Compound> v_cmp;
    vector<ReactionMonitoringTransition> v_rmt_all;
    for (const auto& it : map_mta)
    {
      v_cmp.push_back(it.second.potential_cmp);
      v_rmt_all.insert(v_rmt_all.end(), it.second.potential_rmts.begin(), it.second.potential_rmts.end());
    }

    // convert possible transitions to TargetedExperiment
    TargetedExperiment t_exp;
    t_exp.setCompounds(v_cmp);
    t_exp.setTransitions(v_rmt_all);

    // use MRMAssay methods for filtering
    MRMAssay assay;

    // filter: min/max transitions
    assay.detectingTransitionsCompound(t_exp, min_transitions, max_transitions);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    String extension = out.substr(out.find_last_of(".")+1);

    if (extension == "tsv")
    {
      // validate and write
      OpenMS::TransitionTSVFile::convertTargetedExperimentToTSV(out.c_str(), t_exp);
    }
    else if (extension == "traML")
    {
      // validate
      OpenMS::TransitionTSVFile::validateTargetedExperiment(t_exp);
      // write traML
      TraMLFile traml_out;
      traml_out.store(out, t_exp);
    }
    else if (extension == "pqp")
    {
      //validate 
      OpenMS::TransitionTSVFile::validateTargetedExperiment(t_exp);
      // write pqp
      TransitionPQPFile pqp_out;
      pqp_out.convertTargetedExperimentToPQP(out.c_str(), t_exp);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPAssayGeneratorMetabo tool;
  return tool.main(argc, argv);
}
/// @endcond
