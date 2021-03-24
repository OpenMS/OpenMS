// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>
#include <algorithm>
#include <map>

using namespace OpenMS;

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
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> OpenSWATH pipeline </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_AccurateMassSearch </td>
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

private:
  SiriusAdapterAlgorithm algorithm;

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

    registerDoubleOption_("ambiguity_resolution_mz_tolerance", "<num>", 10.0, "Mz tolerance for the resolution of identification ambiguity over multiple files", false);
    registerStringOption_("ambiguity_resolution_mz_tolerance_unit", "<choice>", "ppm", "Unit of the ambiguity_resolution_mz_tolerance", false, true);
    setValidStrings_("ambiguity_resolution_mz_tolerance_unit", ListUtils::create<String>("ppm,Da"));
    registerDoubleOption_("ambiguity_resolution_rt_tolerance", "<num>", 10.0, "RT tolerance in seconds for the resolution of identification ambiguity over multiple files", false);
    registerDoubleOption_("total_occurrence_filter", "<num>", 0.1, "Filter compound based on total occurrence in analysed samples", false);
    setMinFloat_("total_occurrence_filter", 0.0);
    setMaxFloat_("total_occurrence_filter", 1.0);

    registerDoubleOption_("fragment_annotation_score_threshold", "<num>", 0.80, "Filters annotations based on the explained intensity of the peaks in a spectrum", false);
    setMinFloat_("fragment_annotation_score_threshold", 0.0);
    setMaxFloat_("fragment_annotation_score_threshold", 1.0);

    registerFlag_("decoy_generation", "Decoys will be generated using the fragmentation tree re-rooting approach. This option does only work in combination with the fragment annotation via Sirius.", false);

    registerStringOption_("decoy_generation_method", "<choice>", "original", "Uses different methods for decoy generation. Basis for the method is the fragmentation-tree re-rooting approach ('original'). This approach can be extended by using 'resolve_overlap', which will resolve overlapping fragments of the highest intensity fragments chosen, by adding -CH2 mass to the overlapping fragments. 'Add_shift' will add a -CH2 mass shift to the target fragments and use them as additional decoys if fragmentation-tree re-rooting failed. 'Both' combines the extended methods (resolve_overlap, add_shift).",false);
    setValidStrings_("decoy_generation_method", ListUtils::create<String>("original,resolve_overlap,add_shift,both"));

    registerStringOption_("method", "<choice>", "highest_intensity", "Spectrum with the highest precursor intensity or a consensus spectrum is used for assay library construction (if no fragment annotation is used).",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerFlag_("use_exact_mass", "Use exact mass for precursor and fragment annotations", false);
    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);

    // preprocessing
    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);
    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.1, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("consensus_spectrum_precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds], for consensus spectrum generation (only available without fragment annotation)", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);

    // transition extraction 
    registerIntOption_("min_transitions", "<int>", 3, "Minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "Maximal number of transitions", false);
    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.98, "Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation", false);
    registerDoubleOption_("transition_threshold", "<num>", 5, "Further transitions need at least x% of the maximum intensity (default 5%)", false);
    registerDoubleOption_("min_fragment_mz", "<num>", 0.0, "Minimal m/z of a fragment ion choosen as a transition", false, true);
    registerDoubleOption_("max_fragment_mz", "<num>", 2000.0, "Maximal m/z of a fragment ion choosen as a transition" , false, true);

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
    registerFullParam_(algorithm.getDefaults());
    registerStringOption_("out_workspace_directory", "<directory>", "", "Output directory for SIRIUS workspace", false);
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
    double ar_mz_tol = getDoubleOption_("ambiguity_resolution_mz_tolerance");
    String ar_mz_tol_unit_res = getStringOption_("ambiguity_resolution_mz_tolerance_unit");
    double ar_rt_tol = getDoubleOption_("ambiguity_resolution_rt_tolerance");
    double total_occurrence_filter = getDoubleOption_("total_occurrence_filter");
    double score_threshold = getDoubleOption_("fragment_annotation_score_threshold");
    bool decoy_generation = getFlag_("decoy_generation");
    if (decoy_generation && !use_fragment_annotation)
    {
      decoy_generation = false;
      OPENMS_LOG_WARN << "Warning: Decoy generation was switched off, due to the use of no or an unsupported fragment annotation method." << std::endl;
    }
    bool method_consensus_spectrum = method == "consensus_spectrum" ? true : false;
    bool use_exact_mass = getFlag_("use_exact_mass");
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");

    String decoy_generation_method = getStringOption_("decoy_generation_method");
    bool original = false;
    bool resolve_overlap = false;
    bool add_shift = false;
    if (decoy_generation_method == "original" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting." << std::endl;
      original = true;
    }
    else if (decoy_generation_method == "resolve_overlap" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting and overlap resolution." << std::endl;
      resolve_overlap = true;
    }
    else if (decoy_generation_method == "add_shift" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting and addition of -CH2 mass shift where re-rooting was not possible." << std::endl;
      add_shift = true;
    }
    else if (decoy_generation_method == "both" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting with overlap resolution and addition of -CH2 mass shift where re-rooting was not possible." << std::endl;
      resolve_overlap = true;
      add_shift = true;
    }
    int min_transitions = getIntOption_("min_transitions");
    int max_transitions = getIntOption_("max_transitions");
    double min_fragment_mz = getDoubleOption_("min_fragment_mz");
    double max_fragment_mz = getDoubleOption_("max_fragment_mz");

    double consensus_spectrum_precursor_rt_tolerance = getDoubleOption_("consensus_spectrum_precursor_rt_tolerance");
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

    algorithm.updateExistingParameter(getParam_());

    writeDebug_("Parameters passed to SiriusAdapterAlgorithm", algorithm.getParameters(), 3);

    // SIRIUS workspace (currently needed for fragmentation trees)
    String sirius_workspace_directory = getStringOption_("out_workspace_directory");

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
      FeautreMapping::FeatureMappingInfo fm_info;
      //vector<FeatureMap> v_fp; // copy FeatureMap via push_back
      //KDTreeFeatureMaps fp_map_kd; // reference to *basefeature in vector<FeatureMap>
      FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in vector<FeatureMap>
      algorithm.preprocessingSirius(id[file_counter],
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

      vector< MetaboTargetedAssay::CompoundTargetDecoyPair > v_cmp_spec;

      if (use_fragment_annotation && executable.empty())
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "SIRIUS executable was not found.");
      }
      else if (use_fragment_annotation && !executable.empty())
      {
        // make temporary files
        SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects sirius_tmp(debug_level_);

        // write msfile and store the compound information in CompoundInfo Object
        vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
        SiriusMSFile::store(spectra,
            sirius_tmp.getTmpMsFile(),
            feature_mapping,
            algorithm.isFeatureOnly(),
            algorithm.getIsotopePatternIterations(),
            algorithm.isNoMasstraceInfoIsotopePattern(),
            v_cmpinfo);

        algorithm.logFeatureSpectraNumber(id[file_counter],
            feature_mapping,
            spectra);

        // calls SIRIUS and returns vector of paths to sirius folder structure
        std::vector<String> subdirs;
        String out_csifingerid;
        subdirs = algorithm.callSiriusQProcess(sirius_tmp.getTmpMsFile(),
            sirius_tmp.getTmpOutDir(),
            executable,
            out_csifingerid,
            decoy_generation);
  
        OPENMS_LOG_DEBUG << subdirs.size() << " spectra were annotated using SIRIUS." << std::endl;

        // sort vector path list
        SiriusAdapterAlgorithm::sortSiriusWorkspacePathsByScanIndex(subdirs);

        // extract Sirius/Passatutto FragmentAnnotation and DecoyAnnotation from subdirs
        // and resolve ambiguous identifications in one file based on the native_id_ids and the SIRIUS IsotopeTree_Score
        vector <SiriusFragmentAnnotation::SiriusTargetDecoySpectra> annotated_spectra = SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations(subdirs, score_threshold, use_exact_mass);

        // combine compound information (SiriusMSFile) with annotated spectra (SiriusFragmentAnnotation)
        v_cmp_spec = MetaboTargetedAssay::pairCompoundWithAnnotatedSpectra(v_cmpinfo, annotated_spectra);

        // should the sirius workspace be retained
        if (!sirius_workspace_directory.empty())
        {
          // convert path to absolute path
          QDir sw_dir(sirius_workspace_directory.toQString());
          sirius_workspace_directory = String(sw_dir.absolutePath());

          // move tmp folder to new location
          bool copy_status = File::copyDirRecursively(sirius_tmp.getTmpDir().toQString(), sirius_workspace_directory.toQString());
          if (copy_status)
          {
            OPENMS_LOG_INFO << "Sirius Workspace was successfully copied to " << sirius_workspace_directory << std::endl;
          }
          else
          {
            OPENMS_LOG_INFO << "Sirius Workspace could not be copied to " << sirius_workspace_directory << ". Please run AssayGeneratorMetabo with debug >= 2." << std::endl;
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
                                                                                    min_fragment_mz,
                                                                                    max_fragment_mz,
                                                                                    use_exact_mass,
                                                                                    exclude_ms2_precursor,
                                                                                    file_counter);

      }
      else // use heuristics
      {
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
      }
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

    // sort by highest intensity - filter:  min/max transitions (targets), filter: max transitions (decoys)
    // e.g. if only one decoy fragment is available it will not be filtered out!
    assay.filterMinMaxTransitionsCompound(t_exp, min_transitions, max_transitions);

    // remove decoys which do not have a respective target after min/max transition filtering
    // based on the TransitionGroupID (similar for targets "0_Acephate_[M+H]+_0" and decoys "0_Acephate_decoy_[M+H]+_0")
    if (use_fragment_annotation && decoy_generation)
    {
      assay.filterUnreferencedDecoysCompound(t_exp);
    }

    // resolve overlapping target and decoy masses
    // after selection of decoy masses based on highest intensity (arbitrary, since passatutto uses
    // the intensities based on the previous fragmentation tree), overlapping masses between targets
    // and decoys of one respective metabolite_adduct combination can be resolved by adding a CH2 mass
    if (use_fragment_annotation && decoy_generation && !original)
    {
      const double chtwo_mass = EmpiricalFormula("CH2").getMonoWeight();
      vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> mappings = MetaboTargetedTargetDecoy::constructTargetDecoyMassMapping(t_exp);
      if (resolve_overlap)
      {
        MetaboTargetedTargetDecoy::resolveOverlappingTargetDecoyMassesByIndividualMassShift(t_exp, mappings, chtwo_mass);
      }
      if (add_shift)
      {
        MetaboTargetedTargetDecoy::generateMissingDecoysByMassShift(t_exp, mappings, chtwo_mass);
      }
    }

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
