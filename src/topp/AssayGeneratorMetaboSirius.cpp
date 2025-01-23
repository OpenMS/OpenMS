// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QDir>
#include <QtCore/QDirIterator>
#include <QtCore/QString>
#include <algorithm>
#include <map>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page TOPP_AssayGeneratorMetaboSirius AssayGeneratorMetaboSirius

  @brief Generates an assay library from SIRIUS fragmentation trees (Metabolomics)

    <CENTER>
      <table>
          <tr>
              <th ALIGN = "center"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> &rarr; AssayGeneratorMetaboSirius &rarr;</td>
              <th ALIGN = "center"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> OpenSWATH pipeline </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_AccurateMassSearch </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SiriusExport </td>
          </tr>
      </table>
  </CENTER>

  Build an assay library from a SIRIUS project directory using fragmentation trees.

  - Use the SiriusExport TOPP tool with each samples mzML and featureXML files to generate a SIRIUS input ms file.
  @code
    SiriusExport -in sample.mzML -in_featureinfo sample.featureXML -out_ms sample.ms
  @endcode
  - Run SIRIUS externally with "--no-compression" flag and the formula, passatutto (optional, for decoy generation) and write-summaries tools.
  - This tool was tested with the SIRIUS project directory generated with SIRIUS versions >= 5.5.1 and <= 5.8.6.
  @code
    sirius --input sample.ms --project sirius-project --maxmz 300 --no-compression formula passatutto write-summaries
  @endcode
  - Provide the path to SIRIUS project as input parameter for this tool.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_AssayGeneratorMetaboSirius.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_AssayGeneratorMetaboSirius.html
*/

/// @cond TOPPCLASSES

class TOPPAssayGeneratorMetaboSirius :
  public TOPPBase,
  private TransitionTSVFile
{
public:
  TOPPAssayGeneratorMetaboSirius() :
    TOPPBase("AssayGeneratorMetaboSirius", "Assay library generation from a SIRIUS project directory (Metabolomics)")
    {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<directory>", "", "SIRIUS project directory");

    registerInputFile_("in_compoundinfo", "<file>", "", "Compound info table (.tsv file)");
    setValidFormats_("in_compoundinfo", ListUtils::create<String>("tsv"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv,traML,pqp"));

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


    registerStringOption_("method", "<choice>", "highest_intensity", "Spectrum with the highest precursor intensity or a consensus spectrum is used for assay library construction (if no fragment annotation is used).",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerFlag_("use_exact_mass", "Use exact mass for precursor and fragment annotations", false);
    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);

    // transition extraction
    registerIntOption_("min_transitions", "<int>", 3, "Minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "Maximal number of transitions", false);
    registerDoubleOption_("transition_threshold", "<num>", 5, "Further transitions need at least x% of the maximum intensity (default 5%)", false);
    registerDoubleOption_("min_fragment_mz", "<num>", 0.0, "Minimal m/z of a fragment ion choosen as a transition", false, true);
    registerDoubleOption_("max_fragment_mz", "<num>", 2000.0, "Maximal m/z of a fragment ion choosen as a transition" , false, true);
    
    // decoys
    registerFlag_("decoy_generation", "Decoys will be generated using the fragmentation tree re-rooting approach. This option does only work in combination with the fragment annotation via Sirius.", false);
    registerStringOption_("decoy_generation_method", "<choice>", "original", "Uses different methods for decoy generation. Basis for the method is the fragmentation-tree re-rooting approach ('original'). This approach can be extended by using 'resolve_overlap', which will resolve overlapping target/decoy fragments by adding -CH2 mass to the overlapping decoy fragments. 'generate_missing_decoys' will add a -CH2 mass shift to the target fragments and use them as decoys if fragmentation-tree re-rooting failed. 'Both' combines the extended methods (resolve_overlap, generate_missing_decoys).",false);
    setValidStrings_("decoy_generation_method", ListUtils::create<String>("original,resolve_overlap,generate_missing_decoys,both"));
    registerDoubleOption_("decoy_resolution_mz_tolerance", "<num>", 10.0, "Mz tolerance for the resolution of overlapping m/z values for targets and decoys of one compound.", false);
    registerStringOption_("decoy_resolution_mz_tolerance_unit", "<choice>", "ppm", "Unit of the decoy_resolution_mz_tolerance", false, true);
    setValidStrings_("decoy_resolution_mz_tolerance_unit", ListUtils::create<String>("ppm,Da"));
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------
    String sirius_project_directory = getStringOption_("in");
    String compoundinfo_file = getStringOption_("in_compoundinfo");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    double ar_mz_tol = getDoubleOption_("ambiguity_resolution_mz_tolerance");
    String ar_mz_tol_unit_res = getStringOption_("ambiguity_resolution_mz_tolerance_unit");
    double ar_rt_tol = getDoubleOption_("ambiguity_resolution_rt_tolerance");
    double total_occurrence_filter = getDoubleOption_("total_occurrence_filter");
    double score_threshold = getDoubleOption_("fragment_annotation_score_threshold");
    bool decoy_generation = getFlag_("decoy_generation");
    bool use_exact_mass = getFlag_("use_exact_mass");
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");
    String decoy_generation_method = getStringOption_("decoy_generation_method");
    bool original = false;
    bool resolve_overlap = false;
    bool generate_missing_decoys = false;
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
    else if (decoy_generation_method == "generate_missing_decoys" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting and filling missing decoys by addition of -CH2 mass shift where re-rooting was not possible." << std::endl;
      generate_missing_decoys = true;
    }
    else if (decoy_generation_method == "both" && decoy_generation)
    {
      OPENMS_LOG_INFO << "Decoy method: fragmentation tree re-rooting with overlap resolution and addition of -CH2 mass shift to generate missing decoys where re-rooting was not possible." << std::endl;
      resolve_overlap = true;
      generate_missing_decoys = true;
    }
    double decoy_mz_tol = getDoubleOption_("decoy_resolution_mz_tolerance");
    String decoy_mz_tol_unit_res = getStringOption_("decoy_resolution_mz_tolerance_unit");
    int min_transitions = getIntOption_("min_transitions");
    int max_transitions = getIntOption_("max_transitions");
    double min_fragment_mz = getDoubleOption_("min_fragment_mz");
    double max_fragment_mz = getDoubleOption_("max_fragment_mz");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    //-------------------------------------------------------------
    // Get all subdirectories within the SIRIUS project directory
    //-------------------------------------------------------------
    std::vector<String> subdirs;
    QDirIterator it(sirius_project_directory.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
    while (it.hasNext())
    {
      subdirs.emplace_back(it.next());
    }  
    OPENMS_LOG_DEBUG << subdirs.size() << " spectra were annotated using SIRIUS." << std::endl;
    if (subdirs.empty())
    {
      decoy_generation = false;
      throw OpenMS::Exception::Postcondition(__FILE__,__LINE__, OPENMS_PRETTY_FUNCTION, "SIRIUS project directory is empty.");
    }

    //-------------------------------------------------------------
    // Get CompoundInfo objects from tsv file
    //-------------------------------------------------------------
    vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
    // get number of files from maximum file_index value
    size_t n_files = 0;
    CsvFile csv(compoundinfo_file, '\t');
    size_t row_count = csv.rowCount();
    for (size_t i = 1; i < row_count; ++i) {
      StringList row_data;
      csv.getRow(i, row_data);
      SiriusMSFile::CompoundInfo cmp_info;
      // Convert and assign each field from row_data to cmp_info's attributes
      cmp_info.cmp = row_data[0];
      cmp_info.file_index = stoi(row_data[1]);
      cmp_info.pmass = stod(row_data[2]);
      cmp_info.rt = stod(row_data[4]);
      cmp_info.fmz = stod(row_data[5]);
      cmp_info.fid = row_data[6];
      cmp_info.formula = row_data[7];
      cmp_info.charge = stoi(row_data[8]);
      cmp_info.ionization = row_data[9];
      cmp_info.des = row_data[10];
      cmp_info.source_file = row_data[12];
      cmp_info.m_ids_id = row_data[15];
      // add if "use_known_unknown" flag is set or compound name is not "UNKNOWN"
      if (use_known_unknowns || cmp_info.des != "UNKNOWN") {
          v_cmpinfo.push_back(cmp_info);
      }
      // update n_files with most recent (highest) file_index
      n_files = cmp_info.file_index + 1;
    }

    //--------------------------------------------------------------------------
    // Get list of MetaboTargetedAssay (compound with all possible transitions)
    //--------------------------------------------------------------------------
    //get annotated spectra from SIRIUS project subdirs
    vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra> annotated_spectra =
      SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations(subdirs, score_threshold, use_exact_mass, decoy_generation);

    // combine compound info with annotated spectra
    vector<MetaboTargetedAssay::CompoundTargetDecoyPair> v_cmp_spec;
    v_cmp_spec = MetaboTargetedAssay::pairCompoundWithAnnotatedTDSpectraPairs(v_cmpinfo, annotated_spectra);

    // pair compound info with potential transitions (filtered by min/max, exclude precursor)
    vector<MetaboTargetedAssay> v_mta;
    v_mta = MetaboTargetedAssay::extractMetaboTargetedAssayFragmentAnnotation(v_cmp_spec,
                                                                              transition_threshold,
                                                                              min_fragment_mz,
                                                                              max_fragment_mz,
                                                                              use_exact_mass,
                                                                              exclude_ms2_precursor);

    //--------------------------------------------------------------------------------------------
    // Combine ambigous identifications (derived from consensus features with similar m/z and RT)
    //--------------------------------------------------------------------------------------------
    // build feature maps (matching original raw data files by file_index) and perfom feature linking 
    std::unordered_map< UInt64, vector<MetaboTargetedAssay> > ambiguity_groups = MetaboTargetedAssay::buildAmbiguityGroup(v_mta, ar_mz_tol, ar_rt_tol, ar_mz_tol_unit_res, n_files);

    // resolve identification ambiguity based on highest occurrence and highest intensity
    MetaboTargetedAssay::resolveAmbiguityGroup(ambiguity_groups, total_occurrence_filter, n_files);
    
    //--------------------------------------------------------------------------------------------
    // Merge all transitions in a TargetedExperiment and filter number of transitions
    //--------------------------------------------------------------------------------------------
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

    TargetedExperiment t_exp;
    t_exp.setCompounds(v_cmp);
    t_exp.setTransitions(v_rmt_all);

    // use MRMAssay methods for filtering
    MRMAssay assay;
    // sort by highest intensity - filter:  min/max transitions (targets), filter: max transitions (decoys)
    // e.g. if only one decoy fragment is available it will not be filtered out!
    assay.filterMinMaxTransitionsCompound(t_exp, min_transitions, max_transitions);

    //------------------------------------------------------
    // Decoys
    //------------------------------------------------------
    if (decoy_generation)
    {
      // remove decoys which do not have a respective target after min/max transition filtering
      // based on the TransitionGroupID (similar for targets "0_Acephate_[M+H]+_0" and decoys "0_Acephate_decoy_[M+H]+_0")
      assay.filterUnreferencedDecoysCompound(t_exp);
      // resolve overlapping target and decoy masses
      // after selection of decoy masses based on highest intensity (arbitrary, since passatutto uses
      // the intensities based on the previous fragmentation tree), overlapping masses between targets
      // and decoys of one respective metabolite_adduct combination can be resolved by adding a CH2 mass
      if (!original)
      {
        const double chtwo_mass = EmpiricalFormula("CH2").getMonoWeight();
        vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> mappings = MetaboTargetedTargetDecoy::constructTargetDecoyMassMapping(t_exp);

        if (resolve_overlap)
        {
          MetaboTargetedTargetDecoy::resolveOverlappingTargetDecoyMassesByDecoyMassShift(t_exp, mappings, chtwo_mass, decoy_mz_tol, decoy_mz_tol_unit_res);
        }
        if (generate_missing_decoys)
        {
          MetaboTargetedTargetDecoy::generateMissingDecoysByMassShift(t_exp, mappings, chtwo_mass);
        }
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
  TOPPAssayGeneratorMetaboSirius tool;
  return tool.main(argc, argv);
}
/// @endcond
