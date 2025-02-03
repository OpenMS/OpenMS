// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/RangeUtils.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page TOPP_AssayGeneratorMetabo AssayGeneratorMetabo

  @brief Generates an assay library using DDA data (Metabolomics)

    <CENTER>
      <table>
          <tr>
              <th ALIGN = "center"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> &rarr; AssayGeneratorMetabo &rarr;</td>
              <th ALIGN = "center"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> OpenSWATH pipeline </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_AccurateMassSearch </td>
          </tr>
      </table>
  </CENTER>

  Build an assay library from DDA data (MS and MS/MS) (mzML).
  Please provide a list of features found in the data (featureXML).

  Features can be detected using the FeatureFinderMetabo (FFM) and identifcation information
  can be added using the AccurateMassSearch featureXML output.

  If the FeatureFinderMetabo featureXML does not contain any identifications the "use_known_unknowns" flag is used automatically.

  Internal procedure AssayGeneratorMetabo: \n
  1. Input mzML and featureXML \n
  2. Reannotate precursor mz and intensity \n
  3. Filter feature by number of masstraces \n
  4. Assign precursors to specific feature (FeatureMapping) \n
  5. Extract feature meta information (if possible) \n
  6. Find MS2 spectrum with highest intensity precursor for one feature \n
  7. Annotation is performed either the MS2 with the highest intensity precursor or a consensus spectrum can be used for the transition extraction. \n
  8. Calculate thresholds (maximum and minimum intensity for transition peak) \n
  9. Extract and write transitions (tsv, traml)

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_AssayGeneratorMetabo.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_AssayGeneratorMetabo.html
*/

/// @cond TOPPCLASSES

class TOPPAssayGeneratorMetabo :
  public TOPPBase,
  private TransitionTSVFile
{
public:
  TOPPAssayGeneratorMetabo() :
    TOPPBase("AssayGeneratorMetabo", "Assay library generation from DDA data (Metabolomics)")
    {}

private:
  SiriusExportAlgorithm sirius_export_algorithm;

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file(s)>", StringList(), "MzML input file(s) used for assay library generation");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_featureinfo", "<file(s)>", StringList(), "FeatureXML input file(s) containing identification information (e.g. AccurateMassSearch)");
    setValidFormats_("in_featureinfo", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv,traML,pqp"));

    registerDoubleOption_("ambiguity_resolution_mz_tolerance", "<num>", 10.0, "Mz tolerance for the resolution of identification ambiguity over multiple files", false);
    registerStringOption_("ambiguity_resolution_mz_tolerance_unit", "<choice>", "ppm", "Unit of the ambiguity_resolution_mz_tolerance", false, true);
    setValidStrings_("ambiguity_resolution_mz_tolerance_unit", ListUtils::create<String>("ppm,Da"));
    registerDoubleOption_("ambiguity_resolution_rt_tolerance", "<num>", 10.0, "RT tolerance in seconds for the resolution of identification ambiguity over multiple files", false);
    registerDoubleOption_("total_occurrence_filter", "<num>", 0.1, "Filter compound based on total occurrence in analysed samples", false);
    setMinFloat_("total_occurrence_filter", 0.0);
    setMaxFloat_("total_occurrence_filter", 1.0);

    registerStringOption_("method", "<choice>", "highest_intensity", "Spectrum with the highest precursor intensity or a consensus spectrum is used for assay library construction (if no fragment annotation is used).",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);

    // transition extraction 
    registerIntOption_("min_transitions", "<int>", 3, "Minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "Maximal number of transitions", false);
    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.98, "Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation", false);
    registerDoubleOption_("transition_threshold", "<num>", 5, "Further transitions need at least x% of the maximum intensity (default 5%)", false);
    registerDoubleOption_("min_fragment_mz", "<num>", 0.0, "Minimal m/z of a fragment ion choosen as a transition", false, true);
    registerDoubleOption_("max_fragment_mz", "<num>", 2000.0, "Maximal m/z of a fragment ion choosen as a transition" , false, true);
    

    // precursor
    addEmptyLine_();
    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);
    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.01, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("precursor_consensus_spectrum_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds], for consensus spectrum generation (only available without fragment annotation)", false);

    addEmptyLine_();
    registerFlag_("deisotoping_use_deisotoper", "Use Deisotoper (if no fragment annotation is used)", false);
    registerDoubleOption_("deisotoping_fragment_tolerance", "<num>", 1, "Tolerance used to match isotopic peaks", false);
    registerStringOption_("deisotoping_fragment_unit", "<choice>", "ppm", "Unit of the fragment tolerance", false);
    setValidStrings_("deisotoping_fragment_unit", ListUtils::create<String>("ppm,Da"));
    registerIntOption_("deisotoping_min_charge", "<num>", 1, "The minimum charge considered", false);
    setMinInt_("deisotoping_min_charge", 1);
    registerIntOption_("deisotoping_max_charge", "<num>", 1, "The maximum charge considered", false);
    setMinInt_("deisotoping_max_charge", 1);
    registerIntOption_("deisotoping_min_isopeaks", "<num>", 2, "The minimum number of isotopic peaks (at least 2) required for an isotopic cluster", false);
    setMinInt_("deisotoping_min_isopeaks", 2);
    registerIntOption_("deisotoping_max_isopeaks", "<num>", 3, "The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster", false);
    setMinInt_("deisotoping_max_isopeaks", 3);
    registerFlag_("deisotoping_keep_only_deisotoped", "Only monoisotopic peaks of fragments with isotopic pattern are retained", false);
    registerFlag_("deisotoping_annotate_charge", "Annotate the charge to the peaks", false);

    addEmptyLine_();
    auto defaults = sirius_export_algorithm.getDefaults();
    defaults.remove("isotope_pattern_iterations"); 
    defaults.remove("no_masstrace_info_isotope_pattern"); 

    registerFullParam_(defaults);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    // param AssayGeneratorMetabo
    StringList in = getStringList_("in");
    StringList id = getStringList_("in_featureinfo");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    double ar_mz_tol = getDoubleOption_("ambiguity_resolution_mz_tolerance");
    String ar_mz_tol_unit_res = getStringOption_("ambiguity_resolution_mz_tolerance_unit");
    double ar_rt_tol = getDoubleOption_("ambiguity_resolution_rt_tolerance");
    double total_occurrence_filter = getDoubleOption_("total_occurrence_filter");
    bool method_consensus_spectrum = method == "consensus_spectrum" ? true : false;
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");
    int min_transitions = getIntOption_("min_transitions");
    int max_transitions = getIntOption_("max_transitions");
    double min_fragment_mz = getDoubleOption_("min_fragment_mz");
    double max_fragment_mz = getDoubleOption_("max_fragment_mz");
    double consensus_spectrum_precursor_rt_tolerance = getDoubleOption_("precursor_consensus_spectrum_rt_tolerance");
    double pre_recal_win = getDoubleOption_("precursor_recalibration_window");
    String pre_recal_win_unit = getStringOption_("precursor_recalibration_window_unit");
    bool ppm_recal = pre_recal_win_unit == "ppm" ? true : false;
    double precursor_mz_distance = getDoubleOption_("precursor_mz_distance");
    double cosine_sim_threshold = getDoubleOption_("cosine_similarity_threshold");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    // param deisotoper
    bool use_deisotoper = getFlag_("deisotoping_use_deisotoper");
    double fragment_tolerance = getDoubleOption_("deisotoping_fragment_tolerance");
    String fragment_unit = getStringOption_("deisotoping_fragment_unit");
    bool fragment_unit_ppm = fragment_unit == "ppm" ? true : false;
    int min_charge = getIntOption_("deisotoping_min_charge");
    int max_charge = getIntOption_("deisotoping_max_charge");
    unsigned int min_isopeaks = getIntOption_("deisotoping_min_isopeaks");
    unsigned int max_isopeaks = getIntOption_("deisotoping_max_isopeaks");
    bool keep_only_deisotoped = getFlag_("deisotoping_keep_only_deisotoped");
    bool annotate_charge = getFlag_("deisotoping_annotate_charge");

    writeDebug_("Parameters passed to SiriusExportAlgorithm", sirius_export_algorithm.getParameters(), 3);

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
      PeakMap spectra;
      FileHandler().loadExperiment(in[file_counter], spectra, {FileTypes::MZML});

      // load featurexml
      FeatureMap feature_map;
      FileHandler().loadFeatures(id[file_counter], feature_map, {FileTypes::FEATUREXML});

      // check if featureXML corresponds to mzML
      StringList featurexml_primary_path;
      feature_map.getPrimaryMSRunPath(featurexml_primary_path);

      if (in[file_counter] != featurexml_primary_path[0]) // featureXML should only have one primary path
      {
        OPENMS_LOG_WARN << "Warning: Original paths of the mzML files do not correspond to the featureXML files. Please check and provide the corresponding files." << std::endl;

        OPENMS_LOG_WARN << "Input MzML: " << in[file_counter] << std::endl;

        OPENMS_LOG_WARN << "Input FeatureXML: " << id[file_counter] << std::endl;

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
      FeatureMapping::FeatureMappingInfo fm_info;
      FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in vector<FeatureMap>
      sirius_export_algorithm.preprocessing(id[file_counter],
                                                  spectra,
                                                  fm_info,
                                                  feature_mapping);
    
      // filter known_unkowns based on description (UNKNOWN) (AMS)
      std::map<const BaseFeature*, std::vector<size_t>> feature_ms2_spectra_map = feature_mapping.assignedMS2;
      std::map<const BaseFeature*, std::vector<size_t>> known_features;
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

      // potential transitions of one file
      vector<MetaboTargetedAssay> tmp_mta;
      tmp_mta = MetaboTargetedAssay::extractMetaboTargetedAssay(spectra,
                                                                feature_mapping,
                                                                consensus_spectrum_precursor_rt_tolerance,
                                                                precursor_mz_distance,
                                                                cosine_sim_threshold,
                                                                transition_threshold,
                                                                min_fragment_mz,
                                                                max_fragment_mz,
                                                                method_consensus_spectrum,
                                                                exclude_ms2_precursor,
                                                                file_counter);
      // append potential transitions of one file to vector of all files
      v_mta.insert(v_mta.end(), tmp_mta.begin(), tmp_mta.end());
    } // end iteration over all files

    // group ambiguous identification based on precursor_mz and feature retention time
    // Use featureMap and use FeatureGroupingAlgorithmQT
    std::unordered_map< UInt64, vector<MetaboTargetedAssay> > ambiguity_groups = MetaboTargetedAssay::buildAmbiguityGroup(v_mta, ar_mz_tol, ar_rt_tol, ar_mz_tol_unit_res, in.size());

    // resolve identification ambiguity based on highest occurrence and highest intensity
    MetaboTargetedAssay::resolveAmbiguityGroup(ambiguity_groups, total_occurrence_filter ,in.size());

    // merge possible transitions
    vector<TargetedExperiment::Compound> v_cmp;
    vector<ReactionMonitoringTransition> v_rmt_all;
    for (const auto &it : ambiguity_groups)
    {
      for (const auto &comp_it : it.second)
      {
        v_cmp.push_back(comp_it.potential_cmp);
        v_rmt_all.insert(v_rmt_all.end(), comp_it.potential_rmts.begin(), comp_it.potential_rmts.end());
      }
    }

    // convert possible transitions to TargetedExperiment
    TargetedExperiment t_exp;
    t_exp.setCompounds(v_cmp);
    t_exp.setTransitions(v_rmt_all);

    // use MRMAssay methods for filtering
    MRMAssay assay;

    // sort by highest intensity - filter
    assay.filterMinMaxTransitionsCompound(t_exp, min_transitions, max_transitions);

    // sort TargetedExperiment by name (TransitionID)
    t_exp.sortTransitionsByName();

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
      FileHandler().storeTransitions(out, t_exp, {FileTypes::TRAML});
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
