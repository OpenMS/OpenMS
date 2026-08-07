// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

// the available quantitation methods (instantiated via IsobaricQuantitationMethod::create())
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <memory> // for std::unique_ptr

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/QPXCollectionExport.h>
#include <OpenMS/FORMAT/QPXFile.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IsobaricWorkflow IsobaricWorkflow

@brief Extracts and normalizes isobaric labeling information from an LC-MS/MS experiment.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; IsobaricAnalyzer &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
        </tr>
    </table>
</CENTER>

  The input MSn spectra have to be in centroid mode for the tool to work properly. Use e.g. @ref TOPP_PeakPickerHiRes to perform centroiding of profile data, if necessary.

  This tool currently supports iTRAQ 4-plex and 8-plex, and TMT 6-plex, 10-plex, 11-plex, 16-plex, 18-plex, 32-plex, and 35-plex labeling methods.
  It extracts the isobaric reporter ion intensities from centroided MS2 or MS3 data (MSn), then performs isotope correction and stores the resulting quantitation in a consensus map,
  in which each consensus feature represents one identified PSM together with its reporter ions.
  At least one of @p out, @p out_mzTab, or @p out_qpx must be specified; each output is optional individually.
  The MS level for quantification is chosen automatically per PSM: if MS3 is present, the MS3 product spectrum of the identifying MS2 scan is used (SPS-MS3), otherwise the MS2 scan itself.
  Unlike @ref TOPP_IsobaricAnalyzer, this tool does NOT filter quantification scans by activation method (@p extraction:select_activation is ignored),
  because SPS-MS3 reporter scans are frequently labelled as plain CID rather than HCD; selection is therefore based on the MS-level structure alone.
  For intensity, the closest non-zero m/z signal to the theoretical position is taken as reporter ion abundance.
  The position (RT, m/z) of the consensus centroid is the precursor position in MS1 (from the MS2 spectrum);
  the consensus sub-elements correspond to the theoretical channel m/z (with m/z values of 113-121 Th for iTRAQ and 126-131 Th for TMT, respectively).

  For all labeling techniques, the search radius (@p reporter_mass_shift) should be set as small as possible, to avoid picking up false-positive ions as reporters.
  Usually, Orbitraps deliver precision of about 0.0001 Th at this low mass range. Low intensity reporters might have a slightly higher deviation.
  By default, the mass range is set to ~0.002 Th, which should be sufficient for all instruments (~15 ppm).
  The tool will throw an Exception if you set it below 0.0001 Th (~0.7ppm).
  The tool will also throw an Exception if you set @p reporter_mass_shift > 0.003 Th for TMT-10plex and TMT-11plex, since this could
  lead to ambiguities with neighbouring channels (which are ~0.006 Th apart in most cases).

  For quality control purposes, the tool reports the median distance between the theoretical vs. observed reporter ion peaks in each channel.
  The search radius is fixed to 0.5 Th (regardless of the user defined search radius). This allows to track calibration issues.
  For TMT-10plex, these results are automatically omitted if they could be confused with a neighbouring channel, i.e.
  exceed the tolerance to a neighbouring channel with the same nominal mass (C/N channels).
  If the distance is too large, you might have a m/z calibration problem (see @ref TOPP_InternalCalibration).

  @note If none of the reporter ions can be detected in an MSn scan, a consensus feature will still be generated,
  but the intensities of the overall feature and of all its sub-elements will be zero.
  (If desired, such features can be removed by applying an intensity filter in @ref TOPP_FileFilter.)
  However, if the spectrum is completely empty (no ions whatsoever), no consensus feature will be generated.

  Isotope correction is done using non-negative least squares (NNLS), i.e.:@n
  Minimize ||Ax - b||, subject to x >= 0, where b is the vector of observed reporter intensities (with "contaminating" isotope species),
  A is a correction matrix (as supplied by the manufacturer of the labeling kit) and x is the desired vector of corrected (real) reporter intensities.
  Other software tools solve this problem by using an inverse matrix multiplication, but this can yield entries in x which are negative.
  In a real sample, this solution cannot possibly be true, so usually negative values (= negative reporter intensities) are set to zero.
  However, a negative result usually means that noise was not properly accounted for in the calculation.
  We thus use NNLS to get a non-negative solution, without the need to truncate negative values.
  In the (usual) case that inverse matrix multiplication yields only positive values, our NNLS will give the exact same optimal solution.

  The correction matrices can be found (and changed) in the INI file (parameter @p correction_matrix of the corresponding labeling method).
  However, these matrices for both 4-plex and 8-plex iTRAQ are now stable, and every kit delivered should have the same isotope correction values.
  Thus, there should be no need to change them, but feel free to compare the values in the INI file with your kit's certificate.
  For TMT (6-plex and 10-plex) the values have to be adapted for each kit: Modify the correction matrix according to the data in the product data sheet of your charge:
  <pre>
  Data sheet:
  Mass Tag  Repoter Ion -2      -1      Monoisotopic    +1     +2
  126       126.12776   0.0%    0.0%        100%        5.0%   0.0%
  127N      127.124761  0.0%    0.2%        100%        4.6%   0.0%
  ...
  </pre>
  Corresponding correction matrix:
  <pre>
  [0.0/0.0/5.0/0.0,
  0.0/0.2/4.6/0.0,
  ...
  </pre>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_IsobaricWorkflow.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_IsobaricWorkflow.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIsobaricWorkflow :
  public TOPPBase
{
private:
  std::string ID_RUN_NAME_ = "IsobaricWorkflow_";
  std::map<std::string, std::unique_ptr<IsobaricQuantitationMethod>> quant_methods_;
  std::map<std::string, std::string> quant_method_names_;

  void addMethod_(std::unique_ptr<IsobaricQuantitationMethod> ptr, std::string name)
  {
    std::string internal_name = ptr->getMethodName();
    quant_methods_[internal_name] = std::move(ptr);
    quant_method_names_[internal_name] = name;
  }

public:
  TOPPIsobaricWorkflow() :
    TOPPBase("IsobaricWorkflow", "Calculates isobaric quantitative values for peptides")
  {
    using MT = IsobaricQuantitationMethod::MethodType;
    for (int i = static_cast<int>(MT::UNKNOWN) + 1; i < static_cast<int>(MT::SIZE_OF_METHODTYPE); ++i)
    {
      const MT mt = static_cast<MT>(i);
      addMethod_(IsobaricQuantitationMethod::create(mt), std::string(IsobaricQuantitationMethod::methodDisplayName(mt)));
    }
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // initialize with the first available type
    registerStringOption_("type", "<mode>", quant_methods_.begin()->first, "Isobaric Quantitation method used in the experiment.", false);
    StringList valid_types;
    for (const auto& qm : quant_methods_)
    {
      valid_types.push_back(qm.first);
    }
    setValidStrings_("type", valid_types);

    registerInputFileList_("in", "<file>", {}, "input centroided spectrum files");
    setValidFormats_("in", {"mzML"});
    registerInputFileList_("in_id", "<file>", {}, "corresponding input PSMs");
    setValidFormats_("in_id", {"idXML", "mzId", "idparquet"});
    registerInputFile_("exp_design", "<file>", "", "experimental design file (optional). If not given, the design is assumed to be unfractionated.", false);
    setValidFormats_("exp_design", {"tsv"});
    registerOutputFile_("out", "<file>", "", "Optional output consensusXML file. At least one output must be specified.", false, false);
    setValidFormats_("out", {"consensusXML"});
    registerOutputFile_("out_mzTab", "<file>", "", "Optional output mzTab file with quantitative information. At least one output must be specified.", false, false);
    setValidFormats_("out_mzTab", {"mzTab"});

    registerOutputDir_("out_qpx", "<directory>", "", "Optional output directory for QPX Parquet files (quantms.feature.parquet, quantms.psm.parquet, quantms.pg.parquet). At least one output must be specified.", false, false);
    registerFlag_("calculate_id_purity", "Calculate the purity of the precursor ion based on the MS1 spectrum. Only used for MS3, otherwise it is the same as the quant. precursor purity.");
    registerFlag_("count_sps_matches", "For SPS-MS3: count how many of the co-isolated MS3 precursors (SPS ions) match a b/y fragment ion of the identified peptide. The count is stored as meta value 'sps_matched_ions' on the consensus feature. Off by default.", true);
    registerDoubleOption_("sps_fragment_mass_tolerance", "<tolerance>", 20.0, "Mass tolerance for matching MS3 SPS precursors to theoretical b/y fragment ions of the identified peptide (only used with 'count_sps_matches').", false, true);
    registerStringOption_("sps_fragment_mass_tolerance_unit", "<unit>", "ppm", "Unit of 'sps_fragment_mass_tolerance'.", false, true);
    setValidStrings_("sps_fragment_mass_tolerance_unit", {"ppm", "Da"});
    //registerIntOption_("max_parallel_files", "<num>", 1, "Maximum number of files to load in parallel.", false);
    registerDoubleOption_("psm_score", "<score>", NAN, "The score which should be reached by a peptide hit to be kept.  (use 'NAN' to disable this filter)", false);
    registerDoubleOption_("protein_score", "<score>", NAN, "The score which should be reached by a protein hit to be kept. All proteins are filtered based on their singleton scores irrespective of grouping. Use in combination with 'delete_unreferenced_peptide_hits' to remove affected peptides. (use 'NAN' to disable this filter)", false);
    registerFlag_("delete_unreferenced_peptide_hits", "Peptides not referenced by any protein are deleted in the IDs.");
    // registerFlag_("remove_decoys", "Remove decoys according to the information in the user parameters.");
    registerStringOption_("inference_method", "<option>", "aggregation", "Methods used for protein inference", false);
    setValidStrings_("inference_method", ListUtils::create<std::string>("aggregation,bayesian"));
    registerStringOption_("picked_fdr", "<option>", "false", "Use a picked protein FDR", false, true);
    setValidStrings_("picked_fdr", {"true", "false"});
    registerStringOption_("picked_decoy_string", "<decoy_string>", "", "If using picked protein FDRs, which decoy string was used? Leave blank for auto-detection.", false, true);
    registerStringOption_("picked_decoy_prefix", "<option>", "prefix", "If using picked protein FDRs, was the decoy string a prefix or suffix? Ignored during auto-detection.", false, true);
    setValidStrings_("picked_decoy_prefix", {"prefix", "suffix"});
    registerStringOption_("FDR_type", "<option>", "PSM", "Sub-protein FDR level. PSM, PSM+peptide (best PSM q-value).", false, true);
    setValidStrings_("FDR_type", {"PSM", "PSM+peptide"});
    registerDoubleOption_("proteinFDR", "<threshold>", 1.0, "Protein FDR threshold (0.05=5%).", false);
    setMinFloat_("proteinFDR", 0.0);
    setMaxFloat_("proteinFDR", 1.0);
    registerDoubleOption_("psmFDR", "<threshold>", 1.0, "FDR threshold for sub-protein level (e.g. 0.05=5%)." , false);
    setMinFloat_("psmFDR", 0.0);
    setMaxFloat_("psmFDR", 1.0);
    registerStringOption_("protein_quantification", "<option>", "unique_peptides",
        "Quantify proteins based on:\n"
        "unique_peptides = use peptides mapping to single proteins or a group of indistinguishable proteins"
        "(according to the set of experimentally identified peptides).\n"
        "strictly_unique_peptides = use peptides mapping to a unique single protein only.\n"
        "shared_peptides = use shared peptides only for its best group (by inference score)",
        false, true);
    setValidStrings_("protein_quantification", ListUtils::create<std::string>("unique_peptides,strictly_unique_peptides,shared_peptides"));

    registerSubsection_("extraction", "Parameters for the channel extraction.");
    registerSubsection_("quantification", "Parameters for the peptide quantification.");
    for (const auto& qm : quant_methods_)
    {
      registerSubsection_(qm.second->getMethodName(),std::string("Algorithm parameters for ") + quant_method_names_[qm.second->getMethodName()]);
    }
    Param pq_defaults = PeptideAndProteinQuant().getDefaults();
    pq_defaults.setValue("top:include_all", "true");
    pq_defaults.addTag("top:include_all", "advanced");
    Param bpi_defaults = BasicProteinInferenceAlgorithm().getDefaults();
    // set based on protein_quantification default
    bpi_defaults.remove("annotate_indistinguishable_groups");
    bpi_defaults.remove("greedy_group_resolution");
    Param bayes_defaults = BayesianProteinInferenceAlgorithm().getDefaults();
    
    vector<string> bpi_remove = {"use_ids_outside_features", "annotate_group_probabilities", "user_defined_priors", "update_PSM_probabilities"};
    for (const auto& r : bpi_remove)
    {
      bayes_defaults.remove(r);
    }

    // Why do we only have a const iterator?????
    vector<string> bpi_show = {"psm_probability_cutoff", "top_PSMs"};
    for (auto it = bayes_defaults.begin(); it != bayes_defaults.end(); ++it)
    {
      if (!ListUtils::contains(bpi_show, it.getName()))
      {
        bayes_defaults.addTag(it.getName(), "advanced");
      }
    }

    Param combined;
    combined.insert("ProteinQuantification:", pq_defaults);
    combined.insert("BasicProteinInference:", bpi_defaults);
    combined.insert("BayesianProteinInference:", bayes_defaults);

    registerFullParam_(combined);
  }

  Param getSubsectionDefaults_(const std::string& section) const override
  {
    // any concrete method works to obtain the extractor/quantifier defaults; pick an arbitrary one
    auto temp_quant = IsobaricQuantitationMethod::create(IsobaricQuantitationMethod::MethodType::ITRAQ_4PLEX);
    if (section == "extraction")
    {
      return IsobaricChannelExtractor(temp_quant.get()).getParameters();
    }
    else if (section == "quantification")
    {
      return IsobaricQuantifier(temp_quant.get()).getParameters();
    }
    else
    {
      const auto it = quant_methods_.find(section);
      if (it == quant_methods_.end())
      { // should not happen
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid subsection " + section);
      }
      return it->second->getParameters();
    }
  }

  inline std::tuple<int, int, int> getSpecIdxs_(int pep_idx, const MSExperiment& exp, bool has_ms3)
  {
      int quant_spec_idx = -1;
      int id_spec_idx = pep_idx;
      int ms1_spec_idx = exp.getPrecursorSpectrum(pep_idx);

      if (has_ms3)
      {
        quant_spec_idx = exp.getFirstProductSpectrum(pep_idx);
      }
      else
      {
        quant_spec_idx = pep_idx;
      }

      return std::make_tuple(quant_spec_idx, id_spec_idx, ms1_spec_idx);
  }

  inline std::pair<double, double> getPurities_(
    int quant_spec_idx,
    int id_spec_idx,
    int ms1_spec_idx,
    const MSExperiment& exp,
    bool has_ms3,
    bool max_precursor_isotope_deviation,
    bool calc_id_purity,
    bool interpolate_precursor_purity)
  {
    double quant_purity = -1.0;
    double id_purity = -1.0;
    if (has_ms3)
    {
      // TODO double-check that this is the correct way to compute the purity for MS3. Currently purities are very low.
      std::vector<double> quant_purities = PrecursorPurity::computeSingleScanPrecursorPurities(quant_spec_idx, id_spec_idx, exp, max_precursor_isotope_deviation);
      // average over all precursors
      quant_purity = std::accumulate(quant_purities.begin(), quant_purities.end(), 0.0) / quant_purities.size();
    }
    if (calc_id_purity || !has_ms3)
    {
      double ms1_purity = -1.;
      // An MS2 with no preceding MS1 has no survey scan to measure isolation purity against --
      // an acquisition that starts mid-cycle, or an mzML filtered down to a spectrum subset.
      // MSExperiment::getPrecursorSpectrum() reports that as -1, which is NOT a spectrum index:
      // passing it on reached MSSpectrum::MZBegin() on a spectrum before the start of the
      // experiment and segfaulted inside PrecursorPurity. -1.0 is the value this function
      // already returns for "purity unavailable", and fillConsensusFeature_() only writes the
      // metavalue when it is > 0, so the feature is quantified without a purity annotation
      // instead of the tool dying.
      if (ms1_spec_idx < 0)
      {
        ms1_purity = -1.;
      }
      else if (!interpolate_precursor_purity)
      {
        ms1_purity = PrecursorPurity::computeSingleScanPrecursorPurities(id_spec_idx, ms1_spec_idx, exp, max_precursor_isotope_deviation)[0];
      }
      else
      {
        Size next_ms1_spec = quant_spec_idx;
        do
        {
          next_ms1_spec++;
          if (next_ms1_spec >= exp.size())
          {
            // No more MS1 spectra found, use the original MS1 spectrum
            next_ms1_spec = ms1_spec_idx;
            break;
          }
          else if (exp[next_ms1_spec].getMSLevel() == 1)          
          {
            break;
          }
        } while (next_ms1_spec < exp.size());
        ms1_purity = PrecursorPurity::computeInterpolatedPrecursorPurity(id_spec_idx, ms1_spec_idx, next_ms1_spec, exp, max_precursor_isotope_deviation)[0];
      }

      if (has_ms3)
      {
        id_purity = ms1_purity;
      }
      else
      {
        quant_purity = ms1_purity;
        id_purity = ms1_purity;
      }
    }
    return std::make_pair(quant_purity, id_purity);
  }

  /**
   * @brief Fills an existing ConsensusFeature with all kinds of information of an identified and isobarically quantified peptide.
   *
   * @param[out] cf the ConsensusFeature to fill
   * @param[in] pep information about the PSM (object will be moved)
   * @param[in] exp the MSExperiment to extract information about spectra
   * @param[in] id_spec_idx index of the identifying spectrum
   * @param[in] quant_spec_idx index of the quantifying spectrum
   * @param[in] itys the extracted intensities from the quant. spec.
   * @param[in] quant_method the quantification method used (for channel information), e.g. TMT10plex
   * @param[in] quant_purity purity of the quant. precursor(s) (if available, else -1.)
   * @param[in] id_purity purity of the id precursor
   * @param[in] min_reporter_intensity minimum intensity of a reporter ion to be considered
   * @param[in] file_idx index of the file in the input list
   * @param[in] spec_idx index of the spectrum over all files
   */
  void inline fillConsensusFeature_(ConsensusFeature & cf, PeptideIdentification& pep,
   const MSExperiment& exp, Size id_spec_idx, Size quant_spec_idx, const std::vector<double>& itys,
   const std::unique_ptr<IsobaricQuantitationMethod>& quant_method, double quant_purity, double id_purity,
   double min_reporter_intensity, Size file_idx)
  {
      const auto& quant_spec = exp[quant_spec_idx];
      const auto& id_spec = exp[id_spec_idx];
      //cf.setUniqueId(spec_idx);
      cf.setRT(id_spec.getRT());
      cf.setMZ(id_spec.getPrecursors()[0].getMZ());

      Peak2D channel_value;
      channel_value.setRT(quant_spec.getRT());

      // for each channel of current file
      UInt64 map_index = 0;
      Peak2D::IntensityType overall_intensity = 0.;
      Size col_offset = file_idx * quant_method->getChannelInformation().size();

      for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method->getChannelInformation().begin();
            cl_it != quant_method->getChannelInformation().end();
            ++cl_it)
      {
        // set mz-position of channel
        channel_value.setMZ(cl_it->center);

        // discard contribution of this channel as it is below the required intensity threshold
        if (itys[map_index] < min_reporter_intensity)
        {
          channel_value.setIntensity(0);
        } else {
          channel_value.setIntensity(itys[map_index]);
        }

        overall_intensity += channel_value.getIntensity();

        // add channel to ConsensusFeature
        // TODO for the last param, element_index we probably need a global count of PSMs+channel that we are handling
        //  I hate this useless UniqueIDGeneration
        cf.insert(col_offset + map_index, channel_value, map_index);
        ++map_index;
      } // ! channel_iterator

      // add purity information if we could compute it
      // TODO we should reuse the faster and more efficient quality field of ConsensusFeature
      if (id_purity > 0.0)
      {
        cf.setMetaValue("precursor_purity", id_purity);
      }
      if (quant_purity > 0.0)
      {
        cf.setMetaValue("quant_precursor_purity", quant_purity);
      }

      // embed the id of the scan from which the quantitative information was extracted
      cf.setMetaValue("scan_id", quant_spec.getNativeID());
      // ...as well as additional meta information
      cf.setMetaValue("precursor_intensity", id_spec.getPrecursors()[0].getIntensity());

      cf.setCharge(id_spec.getPrecursors()[0].getCharge());
      cf.setIntensity(overall_intensity);
      pep.setIdentifier(ID_RUN_NAME_);
      // Set id_merge_index to track which input file this peptide identification originated from.
      // This is required for proper mzTab export when multiple files are merged into a single ID run.
      pep.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, file_idx);
      cf.setPeptideIdentifications({std::move(pep)});
  }

  std::string addTimeStamp_(std::string& s) const
  {
    std::array<char, 64> buffer;
    buffer.fill(0);
    time_t rawtime;
    time(&rawtime);
    const auto timeinfo = localtime(&rawtime);
    strftime(buffer.data(), sizeof(buffer), "%d-%m-%Y %H-%M-%S", timeinfo);
    return s + std::string(buffer.data());
  }

  ExitCodes main_(int, const char**) override
  {
    ID_RUN_NAME_ = addTimeStamp_(ID_RUN_NAME_);
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    const std::string out = getStringOption_("out");
    const std::string out_mzTab = getStringOption_("out_mzTab");
    const std::string out_qpx = getOutputDirOption("out_qpx");
    if (out.empty() && out_mzTab.empty() && out_qpx.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                 "out/out_mzTab/out_qpx");
    }

    const std::string exp_design = getStringOption_("exp_design");
    bool bayesian = getStringOption_("inference_method") == "bayesian";
    
    Param pq_param = getParam_().copy("ProteinQuantification:", true);
    writeDebug_("Parameters passed to PeptideAndProteinQuant algorithm", pq_param, 3);

    //-------------------------------------------------------------
    // init quant method and extractor
    //-------------------------------------------------------------
    const auto& quant_method = quant_methods_[getStringOption_("type")];

    // set the parameters for this method
    quant_method->setParameters(getParam_().copy(quant_method->getMethodName() + ":", true));

    bool calc_id_purity = getParam_().getValue("calculate_id_purity").toBool();

    bool count_sps_matches = getParam_().getValue("count_sps_matches").toBool();
    double sps_fragment_mass_tolerance = getDoubleOption_("sps_fragment_mass_tolerance");
    bool sps_fragment_mass_tolerance_unit_ppm = getStringOption_("sps_fragment_mass_tolerance_unit") == "ppm";


    Param extract_param(getParam_().copy("extraction:", true));
    IsobaricChannelExtractor channel_extractor(quant_method.get());
    channel_extractor.setParameters(extract_param);
    double min_reporter_intensity = channel_extractor.getParameters().getValue("min_reporter_intensity");

    // IsobaricWorkflow selects the quantification spectrum per identified PSM purely by MS level
    // (the MS3 product spectrum if MS3 is present, otherwise the identifying MS2 scan). It does NOT
    // filter the quantification scans by activation method: for SPS-MS3 the MS3 reporter scans are
    // frequently labelled as plain CID rather than HCD, so an activation filter would wrongly discard
    // valid reporter scans. If the user explicitly requested a concrete activation method, warn that it
    // is ignored here, so the behaviour is not silently different from IsobaricAnalyzer. See issue #7165.
    {
      const std::string sel_act = channel_extractor.getParameters().getValue("select_activation").toString();
      if (!sel_act.empty() && sel_act != "any" && sel_act != "auto")
      {
        OPENMS_LOG_WARN << "Parameter 'extraction:select_activation' is set to '" << sel_act
                        << "', but IsobaricWorkflow chooses the quantification spectrum automatically by MS level "
                        << "(MS3 if present, otherwise MS2) and does not filter by activation method. "
                        << "This setting will be ignored." << std::endl;
      }
    }

    // TODO since I am mostly using the internal classes IsobaricChannelCorrector and IsobaricNormalizer (if at all),
    //  I should only expose their parameters and only init their objects here.
    IsobaricQuantifier quantifier(quant_method.get());
    Param quant_param(getParam_().copy("quantification:", true));
    quantifier.setParameters(quant_param);

    Matrix<double> correction_matrix = quant_method->getIsotopeCorrectionMatrix();

    // TODO this should go to the purity class or just be a parameter
    bool interpolate_precursor_purity = channel_extractor.getParameters().getValue("purity_interpolation").toBool();
    double max_precursor_isotope_deviation = channel_extractor.getParameters().getValue("precursor_isotope_deviation");

    //const std::string& exp_design = getStringOption_("exp_design");
    IDMergerAlgorithm merger(ID_RUN_NAME_, false);
    ConsensusMap cmap;
    MzMLFile mzml_file;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // iterate over pair of mzML and idXML
    const auto in_mz = getStringList_("in");
    const auto in_id = getStringList_("in_id");
    OPENMS_PRECONDITION(in_mz.size() == in_id.size(), "Number of mzML and idXML files must be equal.");

    /* I tried this but it is the same speed as focussing on the parallelization of the inner loop.
    int max_parallel_files  = std::max((int)getIntOption_("max_parallel_files"), (int)in_mz.size());
    int inner_threads = std::max(1, omp_get_max_threads() / max_parallel_files);
    vector<ConsensusMap> all_cmaps{in_mz.size()};
    #ifdef _OPENMP
    omp_set_max_active_levels(2);
    #pragma omp parallel for num_threads(max_parallel_files)
    #endif
    */

    for (Size i = 0; i < in_mz.size(); ++i)
    {
      //ConsensusMap& cur_cmap = all_cmaps[i];
      ConsensusMap cur_cmap;
      const std::string& mz_file = in_mz[i];
      const std::string& id_file = in_id[i];

      // load mzML
      PeakMap exp;
      mzml_file.load(mz_file, exp);
      std::unordered_map<std::string, Size> ms2scan_to_index;

      bool has_ms3 = false;
      for (Size s = 0; s < exp.size(); ++s)
      {
        if (exp[s].getMSLevel() == 2)
        {
          ms2scan_to_index[exp[s].getNativeID()] = s;
        } else if (exp[s].getMSLevel() == 3)
        {
          has_ms3 = true;
        }
      }

      if(has_ms3)
      {
        OPENMS_LOG_INFO << "Found MS3 spectra. Assuming TMT SPS-MS3 workflow." << std::endl;
      }

      // load idXML
      vector<ProteinIdentification> prot_ids;
      PeptideIdentificationList pep_ids;
      FileHandler().loadIdentifications(id_file, prot_ids, pep_ids,
          {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET}, log_type_);
      // TODO filter by qvalue here?
      double pro_score = getDoubleOption_("protein_score");
      double psm_score = getDoubleOption_("psm_score");

      if (!std::isnan(pro_score)) 
      {
        OPENMS_LOG_INFO << "Filtering by protein score (better than " << pro_score << ")..." << endl;
        IDFilter::filterHitsByScore(prot_ids, pro_score);
      }

      if (!std::isnan(psm_score)) 
      {
        OPENMS_LOG_INFO << "Filtering by PSM score (better than " << psm_score << ")..." << endl;
        IDFilter::filterHitsByScore(pep_ids, psm_score);
      }
      
      merger.insertRuns(std::move(prot_ids), {}); // pep IDs will be stored in the consensus features

      std::vector<ChannelQC> qc;
      qc.resize(quant_method->getNumberOfChannels());
      for (Size i = 0; i < qc.size(); ++i)
      {
        qc[i].mz_deltas = std::vector<double>(pep_ids.size());
      }

      cur_cmap.resize(pep_ids.size());

      channel_extractor.registerChannelsInOutputMap(cmap, mz_file);

      // Collect peptide IDs without corresponding MS3 spectra (thread-safe collection)
      // Note: unassigned_pep_ids is shared across threads; push_back is protected by critical section
      PeptideIdentificationList unassigned_pep_ids;

      #pragma omp parallel for /*num_threads(inner_threads)*/
      for (int64_t pep_idx = 0; pep_idx < static_cast<int64_t>(pep_ids.size()); ++pep_idx)
      {
        auto& pep = pep_ids[pep_idx];
        const auto& spec_ref = pep.getSpectrumReference();
        if (!spec_ref.empty())
        {
          auto ms2spec_it = ms2scan_to_index.find(spec_ref);
          if (ms2spec_it != ms2scan_to_index.end())
          {
            std::vector<std::pair<double,unsigned>> channel_qc(quant_method->getNumberOfChannels(), std::make_pair(std::numeric_limits<double>::quiet_NaN(), 0));

            auto [quant_spec_idx, id_spec_idx, ms1_spec_idx] = getSpecIdxs_(ms2spec_it->second, exp, has_ms3);
            
            // Check if MS3 spectrum is missing for MS2 spectrum
            if (has_ms3 && quant_spec_idx == -1)
            {
              // Store peptide ID with file association for unassigned IDs
              PeptideIdentification unassigned_pep = pep;
              unassigned_pep.setIdentifier(ID_RUN_NAME_);
              unassigned_pep.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, i);
              #pragma omp critical(unassigned_pep_ids_collection)
              {
                OPENMS_LOG_WARN << "MS2 spectrum " << spec_ref << " at index " << ms2spec_it->second
                                << " does not have a corresponding MS3 spectrum. Skipping quantification and adding to unassigned.\n";
                unassigned_pep_ids.push_back(std::move(unassigned_pep));
              }
              continue;
            }
            
            auto [quant_purity, id_purity] = getPurities_(quant_spec_idx, id_spec_idx, ms1_spec_idx, exp, has_ms3, max_precursor_isotope_deviation, calc_id_purity, interpolate_precursor_purity);

            if (has_ms3 && exp[quant_spec_idx].getMSLevel() != 3)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MS3 spectrum expected but not found.",StringUtils::toStr(exp[quant_spec_idx].getMSLevel()));
            }

            std::vector<double> itys = channel_extractor.extractSingleSpec(quant_spec_idx, exp, channel_qc);

            // For SPS-MS3: optionally count how many of the co-isolated MS3 precursors (SPS ions)
            // match a theoretical b/y fragment ion of the identified peptide. This estimates how many
            // of the ions selected for MS3 originate from true fragments of the identified peptide.
            int sps_matched_ions = -1;
            int sps_precursor_count = -1;
            if (count_sps_matches && has_ms3 && !pep.getHits().empty())
            {
              const AASequence& seq = pep.getHits()[0].getSequence();
              const auto& sps_precursors = exp[quant_spec_idx].getPrecursors();
              // SPS fragments of a precursor of charge z are typically singly to (z-1) charged
              const int id_charge = exp[id_spec_idx].getPrecursors().empty() ? 1 : exp[id_spec_idx].getPrecursors()[0].getCharge();
              const int max_fragment_charge = std::max(1, id_charge - 1);
              sps_precursor_count = static_cast<int>(sps_precursors.size());
              sps_matched_ions = static_cast<int>(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments(
                sps_precursors, seq, sps_fragment_mass_tolerance, sps_fragment_mass_tolerance_unit_ppm, max_fragment_charge));
            }

            // TODO if itys are all zero we can actually skip correction and quantification
            // NNLS modifies input, so we need a copy of the correction matrix
            Matrix<double> m = correction_matrix;
            std::vector<double> corrected(itys.size(), 0.);
            NonNegativeLeastSquaresSolver::solve(m, itys, corrected);
            fillConsensusFeature_(cur_cmap[pep_idx], pep, exp, id_spec_idx, quant_spec_idx, corrected, quant_method, quant_purity, id_purity, min_reporter_intensity, i);
            if (sps_matched_ions >= 0)
            {
              cur_cmap[pep_idx].setMetaValue("sps_matched_ions", sps_matched_ions);
              cur_cmap[pep_idx].setMetaValue("sps_precursor_count", sps_precursor_count);
            }
            for (Size i = 0; i < channel_qc.size(); ++i)
            {
              // TODO ProteomeDiscoverer also outputs:
              //  - Reporter S/N but with S/N they mean the intensity corrected for the noise value from the raw thermo Orbitrap files
              //    (and I dont think we read them, see https://github.com/compomics/ThermoRawFileParser/blob/c293d4aa1b04bfd62124ff42c512572427a4316a/Writer/MzMlSpectrumWriter.cs#L1664)
              //  - For SPS-MS3 the number of precursor windows that actually surround a fragment from the identified peptide.
              //    This is available (opt-in) via the 'count_sps_matches' flag and stored as meta value 'sps_matched_ions'.
              qc[i].mz_deltas[pep_idx] = channel_qc[i].first;
              if (channel_qc[i].second > 1)
              {
                #pragma omp atomic
                qc[i].signal_not_unique++;
              }
            }
          }
          else
          {
            // Leaves a default-initialized ConsensusFeature (zero FeatureHandles) in cur_cmap; it is
            // pruned before the merge below. Should not normally happen: an identified spectrum is
            // expected to exist in the mzML.
            OPENMS_LOG_WARN << "Identified spectrum " << spec_ref << " not found in mzML file. Skipping." << std::endl;
          }
        }
      }
      channel_extractor.printStatsWithMissing(qc);

      // Add peptide IDs without MS3 spectra to unassigned IDs
      if (!unassigned_pep_ids.empty())
      {
        OPENMS_LOG_INFO << "Adding " << unassigned_pep_ids.size() 
                        << " peptide identifications without MS3 spectra to unassigned IDs." << std::endl;
        auto& cmap_unassigned = cmap.getUnassignedPeptideIdentifications();
        cmap_unassigned.insert(cmap_unassigned.end(), 
                               std::make_move_iterator(unassigned_pep_ids.begin()), 
                               std::make_move_iterator(unassigned_pep_ids.end()));
      }

      // TODO if we want to support normalization, we either need to replace the quantifier with corrector and normalizer separately
      //  or init the normalizer from the quantifier settings.
      // But honestly, most downstream software can do it better, so I would not bother. Just export to mzTab and do it in R/python.
      //IsobaricNormalizer::normalize(cur_cmap);

      // TODO cleanup, reset?

      // cur_cmap was pre-sized to one ConsensusFeature per peptide ID (resize above) and only the
      // features we could actually quantify were filled in - fillConsensusFeature_ always inserts at
      // least one FeatureHandle per channel. Peptide IDs we could not quantify (an identified MS2 with
      // no corresponding MS3, an identified spectrum missing from the mzML, or one without a spectrum
      // reference) leave their slot as a default-constructed feature with zero FeatureHandles
      // (rt/mz/charge = 0, empty sequence, no map association). Their identifications are preserved in
      // the unassigned peptide IDs, so drop these empty placeholders here instead of merging them into
      // the output map.
      const Size before_prune = cur_cmap.size();
      cur_cmap.erase(remove_if(cur_cmap.begin(), cur_cmap.end(),
                               [](const ConsensusFeature& cf) { return cf.empty(); }),
                     cur_cmap.end());
      if (const Size pruned = before_prune - cur_cmap.size(); pruned > 0)
      {
        OPENMS_LOG_INFO << "Removed " << pruned << " unquantified consensus feature(s) without reporter "
                        << "channels (their peptide identifications are kept in the unassigned IDs)." << std::endl;
      }

      // Always use insert (not move assignment) to preserve column headers registered at line 577.
      // Move assignment would lose column headers if early files have no peptide IDs after filtering.
      cmap.reserve(cmap.size() + cur_cmap.size());
      cmap.insert(cmap.end(), std::make_move_iterator(cur_cmap.begin()), std::make_move_iterator(cur_cmap.end()));

      // TODO If we do the parquet export, we can export the feature file here already. Then, if prot. inference and quant are disabled,
      //  the tool could be run on a single file and distributed over multiple nodes. We could use a parquet partitioned over raw_files
      // If we then implement an inference and quant tool that can read parquet, we can do inference and quant after that parallelized step.
      // Potentially even with pyopenms??
    }

    /*cmap.reserve(std::accumulate(all_cmaps.begin(), all_cmaps.end(), 0, [](Size s, const ConsensusMap& c){return s + c.size();}));
    for (auto& cur_cmap : all_cmaps)
    {
      cmap.insert(cmap.end(), std::make_move_iterator(cur_cmap.begin()), std::make_move_iterator(cur_cmap.end()));
    }*/

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // Annotate output with data processing info
    // Create a custom DataProcessing with filtered parameters
    // (exclude unused quantification method sections to reduce output size)
    DataProcessing dp = getProcessingInfo_(DataProcessing::QUANTITATION);
    
    // Remove parameters for unused quantification methods
    std::string selected_method = quant_method->getMethodName();
    vector<std::string> keys_to_remove;
    
    for (const auto& qm : quant_methods_)
    {
      if (qm.first != selected_method)
      {
        // Collect all parameter keys that start with this unused method name
        vector<std::string> all_keys;
        dp.getKeys(all_keys);
        for (const std::string& key : all_keys)
        {
          if (StringUtils::hasPrefix(key, "parameter: " + qm.first + ":"))
          {
            keys_to_remove.push_back(key);
          }
        }
      }
    }
    for (const std::string& key : keys_to_remove)
    {
      dp.removeMetaValue(key);
    }
    
    addDataProcessing_(cmap, dp);

    // Drop features that carry no reporter signal (zero total intensity). Structurally empty
    // (unquantified) features were already removed per input file before merging above; this remaining
    // filter only affects real, identified features whose reporter channels were all below threshold.
    // TODO expose this as a proper min-intensity filtering setting.
    const auto empty_feat = [](const ConsensusFeature& c){return c.getIntensity() <= 0.;};
    cmap.erase(remove_if(cmap.begin(), cmap.end(), empty_feat), cmap.end());
    cmap.ensureUniqueId();
    std::vector<ProteinIdentification> merged_prot_ids;
    merged_prot_ids.resize(1);
    PeptideIdentificationList _;
    merger.returnResultsAndClear(merged_prot_ids[0], _);
    std::cout << "Merged " << merged_prot_ids[0].getHits().size() << " proteins." << std::endl;
    cmap.setProteinIdentifications(merged_prot_ids);

    // Only fall back to deriving a design from the map when none was supplied. Calling
    // fromConsensusMap() unconditionally logged "No fractions annotated in consensusXML.
    // Assuming unfractionated." once per column header even when -exp_design says otherwise.
    ExperimentalDesign design = (exp_design != "") ? ExperimentalDesignFile::load(exp_design, true)
                                                   : ExperimentalDesign::fromConsensusMap(cmap);

    // Annotate the resolved design back onto the column headers. IsobaricChannelExtractor sets
    // only the channel meta values, so without this the fraction structure is known here and
    // nowhere else -- consumers reading the map cannot tell two fractions of one sample from
    // two independent runs.
    if (const Size unannotated = design.annotateColumnHeaders(cmap); unannotated > 0)
    {
      OPENMS_LOG_WARN << "IsobaricWorkflow: " << unannotated << " of " << cmap.getColumnHeaders().size()
                      << " (file, channel) column(s) have no matching row in the experimental "
                         "design, so they carry no fraction annotation. Downstream fraction-aware "
                         "output cannot group them." << std::endl;
    }

    bool groups = getStringOption_("protein_quantification") != "strictly_unique_peptides";
    bool greedy_group_resolution = getStringOption_("protein_quantification") == "shared_peptides";

    if (!bayesian) {
      BasicProteinInferenceAlgorithm prot_inference;
      Param bpi_param = getParam_().copy("BasicProteinInference:", true);
      bpi_param.setValue("annotate_indistinguishable_groups", groups ? "true" : "false");
      bpi_param.setValue("greedy_group_resolution", greedy_group_resolution ? "true" : "false");
      writeDebug_("Parameters passed to BasicProteinInference algorithm", bpi_param, 3);
      prot_inference.setParameters(bpi_param);
      prot_inference.run(cmap, cmap.getProteinIdentifications()[0], false);
    }
    else {
      BayesianProteinInferenceAlgorithm bayes;
      Param bayes_param = getParam_().copy("BayesianProteinInference:", true);
      // Hardcoded for this usecase
      bayes_param.setValue("update_PSM_probabilities", "false");
      bayes_param.setValue("annotate_group_probabilities", "true");
      bayes_param.setValue("user_defined_priors", "false");
      bayes_param.setValue("use_ids_outside_features", "false");
      writeDebug_("Parameters passed to BayesianProteinInference algorithm", bayes_param, 3);
      bayes.setParameters(bayes_param);
      
      bayes.inferPosteriorProbabilities(cmap, greedy_group_resolution);
      if (!groups) {
        cmap.getProteinIdentifications()[0].getIndistinguishableProteins().clear();
      }
    }

    FalseDiscoveryRate fdr;
    auto& proteins = cmap.getProteinIdentifications()[0];

    if (getStringOption_("picked_fdr") == "true") 
    {
      fdr.applyPickedProteinFDR(proteins, getStringOption_("picked_decoy_string"), getStringOption_("picked_decoy_prefix") == "prefix");
    }
    else 
    {
      fdr.applyBasic(proteins);
    }
    
    if (getStringOption_("FDR_type") == "PSM+peptide")
    { 
      fdr.applyBasicPeptideLevel(cmap, false);
    }
    else 
    {
      fdr.applyBasic(cmap, false);
    }
    

    bool rm_pep = getFlag_("delete_unreferenced_peptide_hits");
    if (rm_pep)
    {
      OPENMS_LOG_INFO << "Removing peptide hits without protein references..." << endl;
    }

    IDFilter::removeDanglingProteinReferences(cmap, rm_pep);
    IDFilter::removeUnreferencedProteins(cmap, true);
    IDFilter::updateProteinGroups(proteins.getIndistinguishableProteins(), proteins.getHits());
    IDFilter::updateProteinGroups(proteins.getProteinGroups(), proteins.getHits());

    const double max_pro_fdr = getDoubleOption_("proteinFDR");
    const double max_psm_fdr = getDoubleOption_("psmFDR");

    // FDR filtering
    if (max_psm_fdr < 1.0) 
    { 
      for (auto& f : cmap) 
      {
        IDFilter::filterHitsByScore(f.getPeptideIdentifications(), max_psm_fdr);
      }
      IDFilter::filterHitsByScore(cmap.getUnassignedPeptideIdentifications(), max_psm_fdr);
    }

    if (max_pro_fdr < 1.0)
    {
      IDFilter::filterHitsByScore(proteins, max_pro_fdr);
      IDFilter::removeDanglingProteinReferences(cmap, rm_pep);
    }

    if (max_psm_fdr < 1.0) 
    { 
      IDFilter::removeUnreferencedProteins(cmap, true); 
    }

    if (max_pro_fdr < 1.0 || max_psm_fdr < 1.0)
    {
      IDFilter::updateProteinGroups(proteins.getIndistinguishableProteins(), proteins.getHits());
      IDFilter::updateProteinGroups(proteins.getProteinGroups(), proteins.getHits());
    }

    if (proteins.getHits().empty())
    {
      throw Exception::MissingInformation(
        __FILE__, 
        __LINE__, 
        OPENMS_PRETTY_FUNCTION,
        "No proteins left after FDR filtering. Please check the log and adjust your settings.");
    }

    if (!greedy_group_resolution && !groups)
    {
      for (auto& f : cmap)
      {
        IDFilter::keepUniquePeptidesPerProtein(f.getPeptideIdentifications());
      }
      IDFilter::keepUniquePeptidesPerProtein(cmap.getUnassignedPeptideIdentifications());
    }


    PeptideAndProteinQuant prot_quantifier;
    prot_quantifier.setParameters(pq_param);
    prot_quantifier.readQuantData(
       cmap,
       design);
    prot_quantifier.quantifyPeptides();
    ProteinIdentification& inferred_proteins = cmap.getProteinIdentifications()[0];
    if (inferred_proteins.getIndistinguishableProteins().empty())
    {
      throw Exception::MissingInformation(
       __FILE__,
       __LINE__,
       OPENMS_PRETTY_FUNCTION,
       "No information on indistinguishable protein groups found.");
    }

    prot_quantifier.quantifyProteins(inferred_proteins);

    auto const & protein_quants = prot_quantifier.getProteinResults();
    if (protein_quants.empty())
    {
     OPENMS_LOG_WARN << "Warning: No proteins were quantified." << endl;
    }

    // Annotate quants to protein(groups) for easier export in mzTab
    // Note: we keep protein groups that have not been quantified
    prot_quantifier.annotateQuantificationsToProteins(
      protein_quants, inferred_proteins, true);

    {
      if (!out_qpx.empty())
      {
        OPENMS_LOG_INFO << "Exporting QPX Parquet files to: " << out_qpx << std::endl;

        // Validate the whole collection before the first write. This tool streams millions of
        // feature and PSM rows before the pg view is even reached, so a design-level refusal
        // there would otherwise land after two very large files are already on disk.
        if (!QPXCollectionExport::requireExportable(cmap, design))
        {
          return CANNOT_WRITE_OUTPUT_FILE; // already logged
        }

        // Whatever the preflight cannot decide up front - a row-level refusal, an I/O error -
        // is undone here: without the commit below, every file already written is removed.
        // This matters most for the streaming views, where a late batch can be refused after
        // gigabytes have been flushed.
        QPXCollectionExport::Transaction qpx_collection(out_qpx);

        // Feature-level export: stream in batches so peak memory stays bounded. For isobaric
        // data there is ~one consensus feature per PSM, so the feature table has millions of
        // rows; the one-shot path builds it all in memory at once and drives large runs into swap.
        // n_threads=0 builds each batch's partitions in parallel across all available cores.
        if (!ConsensusMapArrowExport::exportToParquetStreaming(cmap, out_qpx + "/quantms.feature.parquet",
                                                               /*batch_size=*/1000000,
                                                               ParquetWriteConfig{},
                                                               /*n_threads=*/0))
        {
          OPENMS_LOG_ERROR << "Failed to write features Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        // PSM-level export: collect non-owning pointers to all peptide IDs in the
        // consensus map (assigned per feature, then unassigned). Avoids deep-copying
        // millions of PeptideIdentifications; pointers reference into `cmap`, which
        // outlives the streaming export below.
        std::vector<const PeptideIdentification*> all_pepid_ptrs;
        size_t n_pep = cmap.getUnassignedPeptideIdentifications().size();
        for (const auto& feature : cmap) { n_pep += feature.getPeptideIdentifications().size(); }
        all_pepid_ptrs.reserve(n_pep);
        for (const auto& feature : cmap)
        {
          for (const auto& pepid : feature.getPeptideIdentifications())
          {
            all_pepid_ptrs.push_back(&pepid);
          }
        }
        for (const auto& pepid : cmap.getUnassignedPeptideIdentifications())
        {
          all_pepid_ptrs.push_back(&pepid);
        }

        // Streaming/batched write keeps peak memory bounded for very large PSM counts;
        // n_threads=0 builds each batch's partitions in parallel across all available cores.
        if (!QPXFile::exportToParquetStreaming(cmap.getProteinIdentifications(), all_pepid_ptrs,
                                               out_qpx + "/quantms.psm.parquet",
                                               /*export_all_psms=*/false,
                                               /*batch_size=*/1000000,
                                               ParquetWriteConfig{},
                                               /*n_threads=*/0))
        {
          OPENMS_LOG_ERROR << "Failed to write PSM Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        // Protein group export
        // Pass the design that drove quantification: QPX 1.1 keys the pg view on the set of
        // files aggregated into one quantity, and the design is what defines that grouping and
        // the sample numbering the protein abundances are stored under.
        if (!ProteinGroupArrowExport::exportToParquet(cmap, design, out_qpx + "/quantms.pg.parquet"))
        {
          OPENMS_LOG_ERROR << "Failed to write protein groups Parquet file" << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }

        qpx_collection.commit(); // all three views written
      }
    }

    if (!out.empty())
    {
      FileHandler().storeConsensusFeatures(out, cmap);
    }

    if (!out_mzTab.empty())
    {
      const bool report_unidentified_features(false);
      const bool report_unmapped(true);
      const bool report_subfeatures(false);
      const bool report_unidentified_spectra(false);
      const bool report_not_only_best_psm_per_spectrum(false); 

      MzTabFile().store(out_mzTab, 
                        cmap, 
                        false,
                        report_unidentified_features, 
                        report_unmapped, 
                        report_subfeatures, 
                        report_unidentified_spectra,
                        report_not_only_best_psm_per_spectrum);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPIsobaricWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
