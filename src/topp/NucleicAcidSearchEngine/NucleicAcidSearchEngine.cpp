// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Samuel Wein, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <tuple>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/StatisticFunctions.h> // for "median"
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

// file types
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

// digestion enzymes
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/SYSTEM/File.h>

// ribonucleotides
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include "ModifiedNASequenceGenerator.h"
#include <OpenMS/CHEMISTRY/NASequence.h>

// preprocessing and filtering of spectra
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

// spectra comparison
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

// post-processing of results
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>


#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <regex>

// multithreading
#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_NucleicAcidSearchEngine NucleicAcidSearchEngine

@brief Matches tandem mass spectra to nucleic acid sequences.

Given a FASTA file containing RNA sequences (and optionally decoys) and an mzML file from a nucleic acid mass spec experiment:
- Generate a list of digestion fragments from the FASTA file (based on a specified RNase)
- Search the mzML input for MS2 spectra with parent masses corresponding to any of these sequence fragments
- Match the MS2 spectra to theoretically generated spectra
- Score the resulting matches

Output is in the form of an mzTab-like text file containing the search results.
Optionally, an idXML file suitable for visualizing search results in TOPPView (parameter @p id_out) and a "target coordinates" file for label-free quantification using FeatureFinderMetaboIdent (parameter @p lfq_out) can be generated.

Modified ribonucleotides can either be specified in the FASTA input file if they are expected at a specific site. As globally replacing the unmodified base (as @e fixed modifications), or set as @e variable modifications in the tool options.
Information on available modifications is taken from the Modomics database (http://modomics.genesilico.pl/).
In addition to these "standard" modifications, OpenMS defines "generic" and "ambiguous" ones:
<br>
A generic modification represents a group of modifications that cannot be distinguished by tandem mass spectrometry.
For example, "mA" stands for any methyladenosine (could be "m1A", "m2A", "m6A" or "m8A"), "mmA" for any dimethyladenosine (with two methyl groups on the base), and "mAm" for any 2'-O-dimethyladenosine (with one methyl group each on base and ribose).
There is no technical difference between searching for "mA" or e.g. "m1A", but the generic code better represents that no statement can be made about the position of the methyl group on the base.
<br>
In contrast, an ambiguous modification represents two isobaric modifications (or modification groups) with a methyl group on either the base or the ribose, that could in principle be distinguished based on a-B ions.
For example, "mA?" stands for methyladenosine ("mA", see above) or 2'-O-methyladenosine ("Am").
When using ambiguous modifications in a search, NucleicAcidSearchEngine can optionally try to assign the alternative that generates better a-B ion matches in a spectrum (see parameter @p modifications:resolve_ambiguities).


Exact score ties at @p report:top_hits remain a complete group. Exported ties use
the library's deterministic molecule, charge and adduct value order. A representative selected downstream
is reproducible but does not establish unique modification localization.
When using @p digest, search tolerances and modifications in the output describe
the current run while database and digestion provenance remain from the cache.
LFQ targets use a finite intensity-weighted retention time, falling back to the
median of finite observations when no valid positive intensity weights exist.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_NucleicAcidSearchEngine.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_NucleicAcidSearchEngine.html
*/

class NucleicAcidSearchEngine :
  public TOPPBase
{
  using ConstRibonucleotidePtr = const Ribonucleotide*;

public:
  NucleicAcidSearchEngine() :
    TOPPBase("NucleicAcidSearchEngine", "Annotate nucleic acid identifications to MS/MS spectra."),
    fragment_ion_codes_({"a-B", "a", "b", "c", "d", "w", "x", "y", "z"}),
    resolve_ambiguous_mods_(false)
  {
  }

protected:
  vector<std::string> fragment_ion_codes_;
  map<std::string, std::string> ambiguous_mods_; //< map: specific code -> ambig. code
  bool resolve_ambiguous_mods_;


  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file: spectra");
    setValidFormats_("in", {"mzML",
#ifdef WITH_OPENTIMS
      "d",
#endif
#ifdef WITH_THERMO_RAW
      "raw",
#endif
    });

    registerInputFile_("database", "<file>", "", "Input file: sequence database. Required unless 'digest' is set.", false);
    setValidFormats_("database", ListUtils::create<std::string>("fasta"));

    registerInputFile_("digest", "<file>", "", "Input file: pre-digested sequence database. Can be used instead of 'database'. Sets all 'oligo:...' parameters.", false);
    setValidFormats_("digest", {"idparquet"});

    registerOutputFile_("out", "<file>", "", "Output file: mzTab");
    setValidFormats_("out", ListUtils::create<std::string>("mzTab"));

    registerOutputFile_("id_out", "<file>", "", "Output file: idXML (for visualization in TOPPView)", false);
    setValidFormats_("id_out", ListUtils::create<std::string>("idXML"));

    registerOutputFile_("db_out", "<file>", "", "Output file: native identification bundle (idparquet directory) with the complete search results", false);
    setValidFormats_("db_out", {"idparquet"});

    registerOutputFile_("digest_out", "<file>", "", "Output file: sequence database digest. Ignored if 'digest' input is used.", false);
    setValidFormats_("digest_out", {"idparquet"});

    registerOutputFile_("lfq_out", "<file>", "", "Output file: targets for label-free quantification using FeatureFinderMetaboIdent ('id' input)", false);
    setValidFormats_("lfq_out", vector<std::string>(1, "tsv"));

    registerOutputFile_("theo_ms2_out", "<file>", "", "Output file: theoretical MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("theo_ms2_out", ListUtils::create<std::string>("mzML"));
    registerOutputFile_("exp_ms2_out", "<file>", "", "Output file: experimental MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("exp_ms2_out", ListUtils::create<std::string>("mzML"));

    registerFlag_("decharge_ms2", "Decharge the MS2 spectra for scoring", true);

    registerTOPPSubsection_("precursor", "Precursor (parent ion) options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance (+/- around uncharged precursor mass)", false);

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", ListUtils::create<std::string>("Da,ppm"));

    registerIntOption_("precursor:min_charge", "<num>", -1, "Minimum precursor charge to be considered", false, false);
    registerIntOption_("precursor:max_charge", "<num>", -20, "Maximum precursor charge to be considered", false, false);

    registerFlag_("precursor:include_unknown_charge", "Include MS2 spectra with unknown precursor charge - try to match them in any possible charge between 'min_charge' and 'max_charge', at the risk of a higher error rate", false);

    registerFlag_("precursor:use_avg_mass", "Use average instead of monoisotopic precursor masses (appropriate for low-resolution instruments)", false);

    // Whether to look for precursors with salt adducts
    registerFlag_("precursor:use_adducts", "Consider possible salt adducts (see 'precursor:potential_adducts') when matching precursor masses", false);
    registerStringList_("precursor:potential_adducts", "<list>", ListUtils::create<std::string>("Na:+"), "Adducts considered to explain mass differences. Format: 'Element:Charge(+/-)', i.e. the number of '+' or '-' indicates the charge, e.g. 'Ca:++' indicates +2. Only used if 'precursor:use_adducts' is set.", false, false);

    IntList isotopes = {0, 1, 2, 3, 4};
    registerIntList_("precursor:isotopes", "<list>", isotopes, "Correct for mono-isotopic peak misassignments. E.g.: 1 = precursor may be misassigned to the first isotopic peak. Ignored if 'use_avg_mass' is set.", false, false);

    registerTOPPSubsection_("fragment", "Fragment (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment mass tolerance", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", ListUtils::create<std::string>("Da,ppm"));

    registerStringList_("fragment:ions", "<choice>", fragment_ion_codes_, "Fragment ions to include in theoretical spectra", false);
    setValidStrings_("fragment:ions", fragment_ion_codes_);

    registerOutputFile_("out_preprocessed", "<file>", "", "Optional MS2 spectra after preprocessing", false, true);
    setValidFormats_("out_preprocessed", {"mzML"});
    registerTOPPSubsection_("preprocessing", "Spectrum preprocessing options");
    registerStringOption_("preprocessing:deisotoping", "<method>", "legacy", "Isotope removal: legacy spacing detector or none", false, true);
    setValidStrings_("preprocessing:deisotoping", {"legacy", "none"});
    registerFlag_("preprocessing:filter_window_mower", "Apply WindowMower filter to remove noise peaks by m/z windows", true);
    registerTOPPSubsection_("preprocessing:window_mower", "WindowMower filter parameters");
    registerDoubleOption_("preprocessing:window_mower:windowsize", "<size>", 100.0, "Size of the sliding window along the m/z axis for WindowMower", false, true);
    setMinFloat_("preprocessing:window_mower:windowsize", 1.0);
    registerIntOption_("preprocessing:window_mower:peakcount", "<num>", 50, "Number of peaks that should be kept per window", false, true);
    setMinInt_("preprocessing:window_mower:peakcount", 1);
    registerStringOption_("preprocessing:window_mower:movetype", "<type>", "slide", "Window movement for noise filtering; use jump to reproduce historical NASE preprocessing", false, true);
    setValidStrings_("preprocessing:window_mower:movetype", ListUtils::create<std::string>("slide,jump"));
    registerFlag_("preprocessing:filter_nlargest", "Apply NLargest filter to keep only the top N most intense peaks", true);
    registerTOPPSubsection_("preprocessing:nlargest", "NLargest filter parameters");
    registerIntOption_("preprocessing:nlargest:n", "<num>", 1000, "Number of largest (most intense) peaks to keep per spectrum", false, true);
    setMinInt_("preprocessing:nlargest:n", 1);
    registerFlag_("preprocessing:remove_precursor_peak", "Remove the precursor peak (and isotopic peaks) from the MS2 spectrum", false);
    registerDoubleOption_("preprocessing:precursor_mass_tolerance", "<tolerance>", 1.5, "Tolerance for precursor peak removal (applied to both sides of the peak)", false, true);
    setMinFloat_("preprocessing:precursor_mass_tolerance", 0.0);
    registerStringOption_("preprocessing:precursor_mass_tolerance_unit", "<unit>", "Da", "Unit for precursor peak removal tolerance", false, true);
    setValidStrings_("preprocessing:precursor_mass_tolerance_unit", ListUtils::create<std::string>("Da,ppm"));
    registerIntOption_("preprocessing:precursor_peak_isotopes", "<num>", 2, "Number of isotopic peaks to remove around the precursor (0 = monoisotopic only)", false, true);
    setMinInt_("preprocessing:precursor_peak_isotopes", 0);

    registerTOPPSubsection_("scoring", "Scoring Options");
    registerStringOption_("scoring:method", "<method>", "hyperscore", "Scoring method to use for spectrum matching", false);
    setValidStrings_("scoring:method", ListUtils::create<std::string>("hyperscore,mvh"));
    registerIntOption_("scoring:num_intensity_classes", "<num>", 3, "Number of intensity classes for peak stratification in MVH scoring (only used if scoring:method is 'mvh')", false, true);
    setMinInt_("scoring:num_intensity_classes", 1);
    setMaxInt_("scoring:num_intensity_classes", 7);
    registerDoubleOption_("scoring:tic_fraction", "<fraction>", 0.98, "Fraction of total ion current to retain for MVH scoring (only used if scoring:method is 'mvh')", false, true);
    setMinFloat_("scoring:tic_fraction", 0.5);
    setMaxFloat_("scoring:tic_fraction", 1.0);
    registerStringOption_("scoring:single_peak_credit", "<choice>", "true", "MVH: count at most one match per retained observed peak, using the best mass-accuracy weight. Set false for the previous per-theoretical-ion counting and weighting; all matching annotations are retained in either mode.", false, true);
    setValidStrings_("scoring:single_peak_credit", ListUtils::create<std::string>("true,false"));
    registerFlag_("scoring:use_mass_accuracy", "Apply mass accuracy-based distance weighting to scores (applies to both 'hyperscore' and 'mvh' methods)", true);

    registerDoubleList_("scoring:blacklist_mz", "<values>", DoubleList(), "List of m/z values to exclude from theoretical spectra during scoring (e.g., uninformative c1/d1 ions from unmodified bases). Peaks matching these values within the blacklist tolerance will be removed before scoring.", false, true);
    registerDoubleOption_("scoring:blacklist_tolerance", "<ppm>", 10.0, "Tolerance (in ppm) for matching peaks to blacklisted m/z values", false, true);
    setMinFloat_("scoring:blacklist_tolerance", 0.0);

    registerTOPPSubsection_("modifications", "Modification options");

    // add modified ribos from database
    vector<std::string> all_mods;
    for (const auto& r : *RibonucleotideDB::getInstance())
    {
      if (r->isModified())
      {
        std::string code = r->getCode();
        // commas aren't allowed in parameter string restrictions:
        all_mods.push_back(StringUtils::remove(code, ','));
      }
    }
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<std::string>(""), "Fixed modifications, specified using modomics (https://genesilico.pl/modomics/modifications) terms.", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<std::string>(""), "Variable modifications", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_oligo", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate oligonucleotide", false, false);
    registerFlag_("modifications:resolve_ambiguities", "Attempt to resolve ambiguous modifications (e.g. 'mA?' for 'mA'/'Am') based on a-B ions.\nThis incurs a performance cost because two modifications have to be considered for each case.\nRequires a-B ions to be enabled in parameter 'fragment:ions'.");

    registerTOPPSubsection_("oligo", "Oligonucleotide (digestion) options (ignored if 'digest' input is used)");
    registerIntOption_("oligo:min_size", "<num>", 5, "Minimum size an oligonucleotide must have after digestion to be considered in the search", false);
    registerIntOption_("oligo:max_size", "<num>", 0, "Maximum size an oligonucleotide must have after digestion to be considered in the search, leave at 0 for no limit", false);

    registerIntOption_("oligo:missed_cleavages", "<num>", 1, "Number of missed cleavages", false, false);

    StringList all_enzymes;
    RNaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("oligo:enzyme", "<choice>", "no cleavage", "The enzyme used for RNA digestion", false);
    setValidStrings_("oligo:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Number of top-scoring hits per spectrum to retain; all exact ties at the cutoff are included (0 for all hits)", false, true);
    setMinInt_("report:top_hits", 0);

    registerDoubleOption_("report:require_coverage", "<fraction>", 0.0, "Minimum fraction of internal oligonucleotide positions that must be covered by fragment ions (0.0 = no requirement, 1.0 = full coverage required). Coverage requires at least one fragment of type a, a-B, b, c, d, w, x, y, or z at each position.", false, true);
    setMinFloat_("report:require_coverage", 0.0);
    setMaxFloat_("report:require_coverage", 1.0);

    registerTOPPSubsection_("fdr", "False Discovery Rate options");
    registerStringOption_("fdr:decoy_pattern", "<string>", "", "std::string used as part of the accession to annotate decoy sequences (e.g. 'DECOY_'). Leave empty to skip the FDR/q-value calculation.", false);
    registerDoubleOption_("fdr:cutoff", "<value>", 1.0, "Cut-off for FDR filtering; search hits with higher q-values will be removed", false);
    setMinFloat_("fdr:cutoff", 0.0);
    setMaxFloat_("fdr:cutoff", 1.0);
    registerFlag_("fdr:remove_decoys", "Do not score hits to decoy sequences and remove them when filtering");
  }

  // relevant information about an MS2 precursor ion
  struct PrecursorInfo
  {
    Size scan_index;
    Int charge;
    Size isotope;
    std::optional<AdductInfo> adduct;

    PrecursorInfo(Size scan_index, Int charge, Size isotope, const std::optional<AdductInfo>& adduct = std::nullopt):
        scan_index(scan_index),
        charge(charge),
        isotope(isotope),
        adduct(adduct)
    {
    }
  };

  // Helper function to filter blacklisted m/z values from theoretical spectrum
  void filterBlacklistedIons_(MSSpectrum& theo_spectrum,
                               const vector<double>& blacklist_mz,
                               double tolerance_ppm) const
  {
    if (blacklist_mz.empty()) return;

    vector<Size> indices_to_remove;
    for (Size i = 0; i < theo_spectrum.size(); ++i)
    {
      double mz = theo_spectrum[i].getMZ();
      for (double blacklisted_mz : blacklist_mz)
      {
        double tolerance_da = blacklisted_mz * tolerance_ppm * 1e-6;
        if (abs(mz - blacklisted_mz) <= tolerance_da)
        {
          indices_to_remove.push_back(i);
          break;
        }
      }
    }

    vector<Size> retained;
    for (Size i = 0; i < theo_spectrum.size(); ++i)
    {
      if (!std::binary_search(indices_to_remove.begin(), indices_to_remove.end(), i))
      {
        retained.push_back(i);
      }
    }
    theo_spectrum.select(retained);
  }

  // slimmer structure to store basic hit information
  struct AnnotatedHit
  {
    const RNaseDigestion::DigestedOligo* oligo_ref;
    NASequence sequence;
    double precursor_error_ppm = 0.0; // precursor mass error in ppm
    vector<PeptideHit::PeakAnnotation> annotations; // peak/ion annotations
    const PrecursorInfo* precursor_ref = nullptr; // precursor information
  };

  typedef multimap<double, AnnotatedHit, greater<double>> HitsByScore;

  // query modified residues from database
  set<ConstRibonucleotidePtr> getModifications_(const set<std::string>& mod_names)
  {
    set<ConstRibonucleotidePtr> modifications;
    auto db_ptr = RibonucleotideDB::getInstance();
    std::regex double_digits(R"((\d)(?=\d))");
    for (std::string m : mod_names)
    {
      ConstRibonucleotidePtr mod = 0;
      try
      {
        mod = db_ptr->getRibonucleotide(m);
      }
      catch (Exception::ElementNotFound& /*e*/)
      {
        // commas between numbers were removed - try reinserting them:
        m = std::regex_replace(m, double_digits, "$&,");
        mod = db_ptr->getRibonucleotide(m);
      }
      if (resolve_ambiguous_mods_ && mod->isAmbiguous())
      {
        pair<ConstRibonucleotidePtr, ConstRibonucleotidePtr> alternatives =
          db_ptr->getRibonucleotideAlternatives(m);
        modifications.insert(alternatives.first);
        modifications.insert(alternatives.second);
        // keep track of reverse associations (specific -> ambiguous);
        // constraint: each mod. can only occur in one ambiguity group!
        ambiguous_mods_[alternatives.first->getCode()] = m;
        ambiguous_mods_[alternatives.second->getCode()] = m;
      }
      else
      {
        modifications.insert(mod);
      }
    }
    return modifications;
  }

  /// @brief Count cleavage-sensitive placements already included in a digested candidate.
  Size countCleavageSensitiveMods_(const NASequence& sequence,
                                   const set<std::string>& sensitive_codes) const
  {
    if (sensitive_codes.empty())
    {
      return 0;
    }

    Size count = 0;
    for (const auto& residue : sequence)
    {
      if (sensitive_codes.count(residue.getCode()))
      {
        ++count;
      }
    }
    return count;
  }

  // check for minimum and maximum size
  class HasInvalidLength
  {
    Size min_size_;
    Size max_size_;
  public:
    explicit HasInvalidLength(Size min_size, Size max_size)
      : min_size_(min_size), max_size_(max_size)
    {
    }
    bool operator()(const NASequence& s) const { return (s.size() < min_size_ || s.size() > max_size_); }
  };

  // turn an adduct string (param. "precursor:potential_adducts") into a formula
  // @TODO: adapt "AdductInfo::parseAdductString" and use that
  AdductInfo parseAdduct_(const std::string& adduct)
  {
    StringList parts;
    StringUtils::split(adduct, ':', parts);
    if (parts.size() != 2)
    {
      std::string error = "entry in parameter 'precursor:potential_adducts' does not have two parts separated by ':'";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    error, adduct);
    }

    // determine charge of adduct (by number of '+' or '-')
    Int pos_charge = std::count(parts[1].begin(), parts[1].end(), '+');
    Int neg_charge = std::count(parts[1].begin(), parts[1].end(), '-');
    OPENMS_LOG_DEBUG << ": " << pos_charge - neg_charge << endl;
    if (pos_charge > 0 && neg_charge > 0)
    {
      std::string error = "entry in parameter 'precursor:potential_adducts' mixes positive and negative charges";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    error, adduct);
    }

    EmpiricalFormula ef(parts[0]);
    int charge = (pos_charge > 0) ? pos_charge : -neg_charge;
    return AdductInfo(adduct, ef, charge);
  }

  // spectrum must not contain 0 intensity peaks and must be sorted by m/z
  void deisotopeAndSingleChargeMSSpectrum_(
    MSSpectrum& in,
    Int min_charge,
    Int max_charge,
    double fragment_tolerance,
    bool fragment_unit_ppm,
    bool keep_only_deisotoped = false,
    Size min_isopeaks = 3,
    Size max_isopeaks = 10,
    bool make_single_charged = true)
  {
    if (in.empty()) return;

    MSSpectrum old_spectrum = in;

    // determine charge seeds and extend them
    vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
    vector<Int> features(old_spectrum.size(), -1);
    Int feature_number = 0;

    bool negative_mode = (max_charge < 0);
    Int step = negative_mode ? -1 : 1;

    for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
    {
      double current_mz = old_spectrum[current_peak].getPosition()[0];

      for (Int q = max_charge; abs(q) >= abs(min_charge); q -= step) // important: test charge hypothesis from high to low (in terms of absolute values)
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1) // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / abs(q);
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton) // test for missing peak
            {
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            }
            else
            {
/*
              // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
              Size n_extensions = extensions.size();
              if (n_extensions != 0)
              {
                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
                {
                  if (i < min_isopeaks)
                  {
                    has_min_isopeaks = false;
                  }
                  break;
                }
              }

              // averagine check passed
*/
              extensions.push_back(p);
            }
          }

          if (has_min_isopeaks)
          {
            //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
            mono_isotopic_peak[current_peak] = q;
            for (Size i = 0; i != extensions.size(); ++i)
            {
              features[extensions[i]] = feature_number;
            }
            feature_number++;
          }
        }
      }
    }

    in.clear(false);
    for (Size i = 0; i != old_spectrum.size(); ++i)
    {
      Int z = mono_isotopic_peak[i];
      if (keep_only_deisotoped)
      {
        if (z == 0)
        {
          continue;
        }

        // if already single charged or no decharging selected keep peak as it is
        if (!make_single_charged)
        {
          in.push_back(old_spectrum[i]);
        }
        else // make singly charged
        {
          Peak1D p = old_spectrum[i];
          if (negative_mode) // z < 0 in this case
          {
            z = abs(z);
            p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
          }
          else
          {
            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          }
          in.push_back(p);
        }
      }
      else
      {
        // keep all unassigned peaks
        if (features[i] < 0)
        {
          in.push_back(old_spectrum[i]);
          continue;
        }

        // convert mono-isotopic peak with charge assigned by deisotoping
        if (z != 0)
        {
          if (!make_single_charged)
          {
            in.push_back(old_spectrum[i]);
          }
          else // make singly charged
          {
            Peak1D p = old_spectrum[i];
            if (negative_mode) // z < 0 in this case
            {
              z = abs(z);
              p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
            }
            else
            {
              p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            }
            in.push_back(p);
          }
        }
      }
    }

    in.sortByPosition();
  }


  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra, bool negative_mode, Int min_charge, Int max_charge, bool include_unknown_charge, bool use_window_mower, double window_size, int window_peakcount, const std::string& window_movetype, bool use_nlargest, int nlargest_n, bool remove_precursor, double precursor_mass_tolerance, bool precursor_tolerance_ppm, int precursor_peak_isotopes, const std::string& deisotoping)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    if (use_window_mower)
    {
      Param filter_param = window_mower_filter.getParameters();
      filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
      filter_param.setValue("peakcount", window_peakcount, "The number of peaks that should be kept.");
      filter_param.setValue("movetype", window_movetype, "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
      window_mower_filter.setParameters(filter_param);
    }

    // Note: we expect a higher number for NA than e.g., for peptides
    NLargest nlargest_filter = NLargest(nlargest_n);

    Size n_zero_charge = 0, n_inferred_charge = 0;

#pragma omp parallel reduction(+: n_zero_charge, n_inferred_charge)
    {
#pragma omp for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size();
         ++exp_index)
    {
      MSSpectrum& spec = exp[exp_index];

      // sort by mz
      spec.sortByPosition();

      if (spec.getPrecursors().empty()) continue; // this shouldn't happen
      Int precursor_charge = abs(spec.getPrecursors()[0].getCharge());
      if (precursor_charge == 0) // no charge information
      {
        n_zero_charge++;
        // maybe we are able to infer the charge state:
        if (spec.getPrecursors().size() > 1) // multiplexed PRM experiment
        {
          // all precursors belong to the same parent, but with different charge
          // states; we want to find the precursor with highest charge and infer
          // its charge state:
          map<double, Size> precursors; // precursor: m/z -> index
          for (Size i = 0; i < spec.getPrecursors().size(); ++i)
          {
            precursors[spec.getPrecursors()[i].getMZ()] = i;
          }
          if (precursors.size() >= 2)
          {
            double mz1 = precursors.begin()->first;
            double mz2 = (++precursors.begin())->first;
            double mz_ratio = mz1 / mz2;

            Int step = negative_mode ? -1 : 1;
            Int inferred_charge = 0;
            for (Int charge = max_charge; abs(charge) > abs(min_charge);
                 charge -= step)
            {
              double charge_ratio = (abs(charge) - 1.0) / abs(charge);
              double ratios_ratio = mz_ratio / charge_ratio;
              if ((ratios_ratio > 0.99) && (ratios_ratio < 1.01))
              {
                inferred_charge = charge;
                break;
              }
            }
            if (inferred_charge == 0)
            {
              OPENMS_LOG_ERROR
                << "Error: unable to determine charge state for spectrum '"
                << spec.getNativeID() << "' based on precursor m/z values "
                << mz1 << " and " << mz2 << endl;
            }
            else
            {
              ++n_inferred_charge;
              OPENMS_LOG_DEBUG << "Inferred charge state " << inferred_charge
                               << " for spectrum '" << spec.getNativeID() << "'"
                               << endl;
              // keep only precursor with highest charge, set inferred charge:
              Precursor prec = spec.getPrecursors()[precursors.begin()->second];
              prec.setCharge(abs(inferred_charge));
              spec.setPrecursors(vector<Precursor>(1, prec));
            }
          }
        }
      }

      precursor_charge = abs(spec.getPrecursors()[0].getCharge());

      // Remove in the original charge/mass space, before decharging or filtering.
      if (remove_precursor)
      {
        const double precursor_mz = spec.getPrecursors()[0].getMZ();
        const Int charge = std::max(1, precursor_charge);
        const double tolerance = precursor_tolerance_ppm ?
          precursor_mz * precursor_mass_tolerance * 1e-6 : precursor_mass_tolerance;
        vector<Size> retained;
        for (Size i = 0; i < spec.size(); ++i)
        {
          bool remove = false;
          for (Int isotope = 0; isotope <= precursor_peak_isotopes; ++isotope)
          {
            const double isotope_mz = precursor_mz + isotope * Constants::C13C12_MASSDIFF_U / charge;
            if (abs(spec[i].getMZ() - isotope_mz) <= tolerance)
            {
              remove = true;
              break;
            }
          }
          if (!remove) retained.push_back(i);
        }
        spec.select(retained);
      }

      // deisotope
      Int coef = negative_mode ? -1 : 1;
      // @TODO: what happens here if "precursor_charge" is zero?
      if (deisotoping == "legacy")
      {
        deisotopeAndSingleChargeMSSpectrum_(spec, coef, coef * precursor_charge, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 20, single_charge_spectra);
      }
      // remove noise
      if (use_window_mower)
      {
        window_mower_filter.filterPeakSpectrum(spec);
      }
      if (use_nlargest)
      {
        nlargest_filter.filterPeakSpectrum(spec);
      }

      // sort (nlargest changes order)
      spec.sortByPosition();
    }

    } // parallel region

    if (n_zero_charge)
    {
      OPENMS_LOG_WARN << "Warning: no charge state information available for "
                      << n_zero_charge << " out of " << exp.size()
                      << " spectra." << endl;
      if (n_inferred_charge)
      {
        OPENMS_LOG_INFO << "Inferred charge states for " << n_inferred_charge
                        << " spectra." << endl;
      }
      if (n_zero_charge - n_inferred_charge > 0)
      {
        OPENMS_LOG_INFO
          << "Spectra without charge information will be "
          << (include_unknown_charge ? "included in the processing" : "skipped")
          << " (see parameter 'precursor:include_unknown_charge')" << endl;
      }
    }
  }


  double calculatePrecursorMass_(double mz, Int charge, Int isotope,
                                 double adduct_mass, bool negative_mode)
  {
    // we want to calculate the unadducted (!) precursor mass at neutral charge:
    double mass = mz * charge - adduct_mass;
    // compensate for loss or gain of protons that confer the charge:
    if (negative_mode)
    {
      mass += Constants::PROTON_MASS_U * charge;
    }
    else
    {
      mass -= Constants::PROTON_MASS_U * charge;
    }
    // correct for precursor not being the monoisotopic peak:
    mass -= isotope * Constants::C13C12_MASSDIFF_U;

    return mass;
  }


  void resolveAmbiguousMods_(HitsByScore& hits)
  {
    OPENMS_PRECONDITION(hits.size() > 1, "more than one hit expected");
    // If the current hit is an ambiguity variant of an earlier hit in the same
    // hypothesis, combine both into one hit. For example, these sequences:
    // 1. "AUC[mA]Gp", 2. "AUC[Am]Gp"
    // The result should be: 1. "AUC[mA?]Gp" (note ambiguity code), 2. removed.
    const auto merge_variants = [&](AnnotatedHit& previous_hit, const AnnotatedHit& current_hit)
    {
      NASequence& previous_seq = previous_hit.sequence;
      const NASequence& current_seq = current_hit.sequence;
      if ((current_seq.size() != previous_seq.size()) ||
          (current_seq.getFivePrimeMod() != previous_seq.getFivePrimeMod()) ||
          (current_seq.getThreePrimeMod() != previous_seq.getThreePrimeMod()) ||
          (current_hit.oligo_ref != previous_hit.oligo_ref) ||
          (current_hit.precursor_ref != previous_hit.precursor_ref)) // different hypotheses
      {
        return false;
      }
      bool remove_current = true;
      NASequence replacement; // potential replacement sequence for previous hit
      for (Size i = 0; i < current_seq.size(); ++i)
      {
        if (previous_seq[i]->getCode() == current_seq[i]->getCode()) continue;
        if (const auto pos_current = ambiguous_mods_.find(current_seq[i]->getCode()); pos_current == ambiguous_mods_.end())
        {
          // difference is not due to an ambiguous mod. - don't combine hits:
          remove_current = false;
          break;
        }
        else
        {
          // is this ribonucleotide in the previous hit already ambiguous?
          const std::string& ambig_code = pos_current->second;
          if (previous_seq[i]->getCode() == ambig_code) continue;
          // if not, should we replace it with an ambiguous mod.?
          if (const auto pos_previous = ambiguous_mods_.find(previous_seq[i]->getCode()); 
              (pos_previous == ambiguous_mods_.end()) ||
              (pos_previous->second != ambig_code)) // mods don't match
          {
            remove_current = false;
            break;
          }
          if (replacement.empty()) replacement = previous_seq;
          replacement[i] = RibonucleotideDB::getInstance()->
            getRibonucleotide(ambig_code);
        }
      }
      if (remove_current) // current hit is redundant -> remove it
      {
        if (!replacement.empty()) previous_seq = replacement;
        return true;
      }
      return false;
    };

    // Sequence sorting can interleave different precursors and unrelated tied
    // sequences. Keep all earlier representatives for each exact-score hypothesis.
    map<const RNaseDigestion::DigestedOligo*,
        map<const PrecursorInfo*, vector<HitsByScore::iterator>>> representatives;
    double score = hits.begin()->first;
    for (auto hit_it = hits.begin(); hit_it != hits.end(); /* no ++ here! */)
    {
      if (hit_it->first != score)
      {
        representatives.clear();
        score = hit_it->first;
      }
      auto& earlier_hits = representatives[hit_it->second.oligo_ref][hit_it->second.precursor_ref];
      const bool merged = any_of(earlier_hits.begin(), earlier_hits.end(), [&](const auto& previous_it)
      {
        return merge_variants(previous_it->second, hit_it->second);
      });
      if (merged)
      {
        hit_it = hits.erase(hit_it);
      }
      else
      {
        earlier_hits.push_back(hit_it);
        ++hit_it;
      }
    }
  }


  /// @brief Sort tied hits by molecule, precursor, and adduct values before registration.
  void sortHitsByValue_(HitsByScore& hits) const
  {
    vector<pair<double, AnnotatedHit>> ordered_hits(hits.begin(), hits.end());
    auto hit_key = [](const AnnotatedHit& hit)
    {
      const auto& precursor = *hit.precursor_ref;
      const std::string adduct = precursor.adduct ? precursor.adduct->getName() : std::string();
      return make_tuple(hit.sequence.toString(), hit.oligo_ref->sequence.toString(),
                        precursor.charge, precursor.isotope, adduct,
                        precursor.adduct ? precursor.adduct->getEmpiricalFormula().toString() : std::string(),
                        precursor.adduct ? precursor.adduct->getCharge() : 0,
                        precursor.adduct ? precursor.adduct->getMolMultiplier() : UInt(0),
                        hit.precursor_error_ppm);
    };
    sort(ordered_hits.begin(), ordered_hits.end(), [&](const auto& lhs, const auto& rhs)
    {
      if (lhs.first != rhs.first)
      {
        return lhs.first > rhs.first;
      }
      return hit_key(lhs.second) < hit_key(rhs.second);
    });
    hits.clear();
    hits.insert(ordered_hits.begin(), ordered_hits.end());
  }

  /// @brief Retain complete score groups at the report cutoff, independent of arrival order.
  void retainTopHit_(HitsByScore& hits, double score, AnnotatedHit hit, Size top_hits) const
  {
    if ((top_hits == 0) || (hits.size() < top_hits) || (score >= hits.rbegin()->first))
    {
      hits.emplace(score, std::move(hit));
      if (top_hits > 0)
      {
        auto cutoff = hits.begin();
        std::advance(cutoff, std::min<Size>(top_hits, hits.size()) - 1);
        hits.erase(hits.upper_bound(cutoff->first), hits.end());
      }
    }
  }

  void postProcessHits_(const PeakMap& exp, vector<HitsByScore>& annotated_hits, IdentificationData& id_data, bool negative_mode)
  {
    auto& run = id_data.getRun("NASE");
    // Scoring is parallel. Construct owning records in spectrum order with one writer.
    for (Size scan_index = 0; scan_index < annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) continue;
      const auto& spectrum = exp[scan_index];
      IdentificationData::Observation observation;
      observation.data_id = spectrum.getNativeID();
      observation.rt = spectrum.getRT();
      observation.mz = spectrum.getPrecursors()[0].getMZ();
      observation.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
      observation.setMetaValue("precursor_intensity", spectrum.getPrecursors()[0].getIntensity());
      auto query = run.addIdentification(run.getSourceId(0), observation);
      sortHitsByValue_(annotated_hits[scan_index]);
      if (resolve_ambiguous_mods_ && annotated_hits[scan_index].size() > 1)
      {
        resolveAmbiguousMods_(annotated_hits[scan_index]);
        // Merged ambiguity codes can change the value order of the surviving hits.
        sortHitsByValue_(annotated_hits[scan_index]);
      }
      for (const auto& [score, hit] : annotated_hits[scan_index])
      {
        IdentificationData::MatchData match;
        match.encoding = IdentificationData::Encoding::NA_SEQUENCE;
        match.representation = hit.sequence.toString();
        match.sequence_evidence = hit.oligo_ref->sequence_evidence;
        match.target_decoy = hit.oligo_ref->target_decoy;
        match.charge = negative_mode ? -std::abs(hit.precursor_ref->charge) : hit.precursor_ref->charge;
        match.peak_annotations = hit.annotations;
        if (hit.precursor_ref->adduct)
        {
          const auto& salt = *hit.precursor_ref->adduct;
          auto formula = salt.getEmpiricalFormula();
          formula += EmpiricalFormula("H") * (match.charge - salt.getCharge());
          const auto& adduct = match.details.emplace().adduct = AdductInfo(salt.getName(), formula, match.charge, salt.getMolMultiplier());
          match.calculated_mz = adduct->getMZ(hit.sequence.getMonoWeight());
        }
        else
          match.calculated_mz = hit.sequence.getMonoWeight(NASequence::Full, match.charge) / std::abs(match.charge);
        match.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, hit.precursor_error_ppm);
        match.setMetaValue("isotope_offset", hit.precursor_ref->isotope);
        run.addMatch(query, match, {score});
      }
    }
    id_data.validate();
  }


  /// Check if all internal positions of an oligonucleotide have fragment coverage
  double calculateFragmentCoverage_(const NASequence& oligo,
                                     const vector<PeptideHit::PeakAnnotation>& annotations)
  {
    Size length = oligo.size();
    if (length <= 1) return 1.0; // No internal positions to check - full coverage

    // Track which positions have coverage (positions 1 through length-1)
    set<Size> covered_positions;

    for (const auto& annotation : annotations)
    {
      std::string ann = annotation.annotation;
      if (ann.empty())
      {
        OPENMS_LOG_WARN << "Empty fragment annotation found in coverage check" << endl;
        continue;
      }

      // Check for a-B ions (e.g., "a1-B", "a2-B")
      if (ann.find("-B") != std::string::npos)
      {
        // Extract position from "a1-B" -> "1"
        std::string prefix = ann.substr(0, ann.find('-'));
        if (prefix.starts_with("a"))
        {
          std::string pos_str = prefix.substr(1); // Remove 'a'
          try
          {
            Size pos = static_cast<Size>(StringUtils::toInt64(pos_str));
            if (pos >= 1 && pos < length)
            {
              covered_positions.insert(pos);
            }
          }
          catch (Exception::ConversionError&)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Malformed a-B fragment annotation: '" + ann + "' - expected format 'aN-B' where N is a position number",
              ann);
          }
        }
      }
      else
      {
        // Regular ions: extract ion type and position
        // Format: "a1", "b2", "w9", "y3", etc.
        char ion_type = ann[0];
        std::string ion_type_str(1, ion_type);

        // Check if this is a valid fragment ion type (using class member)
        if (std::find(fragment_ion_codes_.begin(), fragment_ion_codes_.end(),
                      ion_type_str) != fragment_ion_codes_.end())
        {
          // Extract numeric position
          std::string pos_str = ann.substr(1);
          // Handle annotations like "y3+" by removing non-digits
          pos_str.erase(std::remove(pos_str.begin(), pos_str.end(), '+'), pos_str.end());
          pos_str.erase(std::remove(pos_str.begin(), pos_str.end(), '-'), pos_str.end());

          try
          {
            Size pos = static_cast<Size>(StringUtils::toInt64(pos_str));

            // w/x/y/z ions count from 3' end, need to convert to 5' end position
            // For a 10-mer: w9 = position 1, w8 = position 2, ..., w1 = position 9
            Size actual_pos = pos;
            if (ion_type == 'w' || ion_type == 'x' || ion_type == 'y' || ion_type == 'z')
            {
              actual_pos = length - pos;
            }

            if (actual_pos >= 1 && actual_pos < length)
            {
              covered_positions.insert(actual_pos);
            }
          }
          catch (Exception::ConversionError&)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Malformed fragment annotation: '" + ann + "' - expected format like 'a1', 'b2', 'y3', etc.",
              ann);
          }
        }
      }
    }

    // Calculate coverage fraction: number of covered positions / total internal positions
    Size total_internal_positions = length - 1;
    double coverage = static_cast<double>(covered_positions.size()) / static_cast<double>(total_internal_positions);

    return coverage;
  }


  void calculateAndFilterFDR_(IdentificationData& id_data, bool only_top_hits)
  {
    const auto score_ref = id_data.getRun("NASE").getScoreDefinitions().at(0);
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getDefaults();
    fdr_params.setValue("use_all_hits", only_top_hits ? "false" : "true");
    bool remove_decoys = getFlag_("fdr:remove_decoys");
    fdr_params.setValue("add_decoy_peptides", remove_decoys ? "false" : "true");
    fdr.setParameters(fdr_params);
    auto fdr_ref = fdr.applyToObservationMatches(id_data, score_ref);
    double fdr_cutoff = getDoubleOption_("fdr:cutoff");
    if (remove_decoys) // remove references to decoys from shared oligos
    {
      IDFilter::removeDecoys(id_data);
    }
    if (fdr_cutoff < 1.0)
    {
      IDFilter::filterObservationMatchesByScore(id_data, fdr_ref, fdr_cutoff);
      OPENMS_LOG_INFO << "Search hits after FDR filtering: " << id_data.getRun("NASE").getNumberOfMatches()
                      << "\nIdentified spectra after FDR filtering: " << id_data.getRun("NASE").getNumberOfIdentifications() << endl;
    }
  }


  void generateLFQInput_(IdentificationData& id_data, const std::string& out_file)
  {
    using AdductedOligo = pair<NASequence, std::string>;
    using PrecursorPair = pair<double, double>; // precursor intensity, RT
    // mapping: charge -> list of precursors
    using PrecursorsByCharge = map<Int, vector<PrecursorPair>>;
    map<AdductedOligo, PrecursorsByCharge> rt_info;
    for (const auto& run : id_data.getRuns())
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
          {
            if (! query.rt) continue;
            auto seq = NASequence::fromString(match.representation);
            const auto& adduct = match.details.value_or_default().adduct;
            auto key = make_pair(seq, adduct ? adduct->getName() : std::string {});
            double intensity = query.getMetaValue("precursor_intensity");
            rt_info[key][match.charge].push_back({intensity, *query.rt});
          }

    SVOutStream tsv(out_file);
    tsv.modifyStrings(false);
    tsv << "CompoundName" << "SumFormula" << "Mass" << "Charge"
        << "RetentionTime" << "RetentionTimeRange" << "IsoDistribution" << endl;
    for (const auto& entry : rt_info)
    {
      std::string name = entry.first.first.toString();
      EmpiricalFormula ef = entry.first.first.getFormula();
      const std::string& adduct = entry.first.second;
      if (! adduct.empty())
      {
        name += "+[" + adduct + "]";
        ef += parseAdduct_(adduct).getEmpiricalFormula();
      }
      // @TODO: use charge-specific RTs?
      vector<Int> charges;
      vector<double> rts;
      for (const auto& charge_pair : entry.second)
      {
        // use intensity-weighted mean of precursor RTs as "apex" RT:
        double weighted_rt = 0.0, total_weight = 0.0;
        for (const auto& rt_pair : charge_pair.second)
        {
          if (std::isfinite(rt_pair.first) && rt_pair.first > 0.0 && std::isfinite(rt_pair.second))
          {
            weighted_rt += rt_pair.first * rt_pair.second;
            total_weight += rt_pair.first;
          }
        }
        if (total_weight > 0.0 && std::isfinite(total_weight) && std::isfinite(weighted_rt))
        {
          rts.push_back(weighted_rt / total_weight);
          charges.push_back(charge_pair.first);
        }
        else
        {
          vector<double> observed_rts;
          for (const auto& rt_pair : charge_pair.second)
          {
            if (std::isfinite(rt_pair.second))
            {
              observed_rts.push_back(rt_pair.second);
            }
          }
          if (!observed_rts.empty())
          {
            rts.push_back(Math::median(observed_rts.begin(), observed_rts.end(), false));
            charges.push_back(charge_pair.first);
          }
        }
      }
      if (rts.empty())
      {
        OPENMS_LOG_WARN << "Skipping LFQ target with no finite retention time: " << name << endl;
        continue;
      }
      tsv << name << ef << 0 << ListUtils::concatenate(charges, ",");
      // overall target RT is median over all charge states:
      tsv << Math::median(rts.begin(), rts.end(), false) << 0 << 0 << endl;
    }
  }


  ExitCodes main_(int, const char**) override
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    IdentificationData id_data;
    auto& result_run = id_data.addRun("NASE", IdentificationData::MoleculeKind::OLIGONUCLEOTIDE);
    vector<RNaseDigestion::DigestedOligo> digest;

    // load parameters and check validity:
    std::string in_mzml = getStringOption_("in");
    std::string in_db = getStringOption_("database");
    std::string in_digest = getStringOption_("digest");

    if (in_db.empty() && in_digest.empty())
    {
      OPENMS_LOG_ERROR << "Error: parameter 'database' or 'digest' must be set"
                       << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (!in_db.empty() && !in_digest.empty())
    {
      OPENMS_LOG_WARN
        << "Warning: both 'database' and 'digest' are set; ignoring 'database'"
        << endl;
    }

    std::string out = getStringOption_("out");
    std::string id_out = getStringOption_("id_out");
    std::string db_out = getStringOption_("db_out");
    std::string lfq_out = getStringOption_("lfq_out");
    std::string theo_ms2_out = getStringOption_("theo_ms2_out");
    std::string exp_ms2_out = getStringOption_("exp_ms2_out");
    bool use_avg_mass = getFlag_("precursor:use_avg_mass");
    Int min_charge = getIntOption_("precursor:min_charge");
    Int max_charge = getIntOption_("precursor:max_charge");

    // @TODO: allow zero to mean "any charge state in the data"?
    if ((min_charge == 0) || (max_charge == 0))
    {
      OPENMS_LOG_ERROR << "Error: invalid charge state 0" << endl;
      return ILLEGAL_PARAMETERS;
    }
    // charges can be positive or negative, depending on data acquisition mode:
    if (((min_charge < 0) && (max_charge > 0)) ||
        ((min_charge > 0) && (max_charge < 0)))
    {
      OPENMS_LOG_ERROR << "Error: mixing positive and negative charges is not allowed"
                       << endl;
      return ILLEGAL_PARAMETERS;
    }
    // min./max. are based on absolute value:
    if (abs(max_charge) < abs(min_charge)) swap(min_charge, max_charge);
    bool negative_mode = (max_charge < 0);
    Int charge_step = negative_mode ? -1 : 1;

    struct SearchConfiguration
    {
      set<Int> charges;
      double precursor_mass_tolerance = 0, fragment_mass_tolerance = 0;
      bool precursor_tolerance_ppm = false, fragment_tolerance_ppm = false;
      Size min_length = 0, max_length = 0, missed_cleavages = 0;
      set<std::string> fixed_mods, variable_mods;
      std::string database;
      const DigestionEnzyme* digestion_enzyme = nullptr;
    } search_param;
    for (Int charge = min_charge; abs(charge) <= abs(max_charge);
         charge += charge_step)
    {
      search_param.charges.insert(charge);
    }
    search_param.precursor_mass_tolerance =
      getDoubleOption_("precursor:mass_tolerance");
    search_param.precursor_tolerance_ppm =
      (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    search_param.fragment_mass_tolerance =
      getDoubleOption_("fragment:mass_tolerance");
    search_param.fragment_tolerance_ppm =
      (getStringOption_("fragment:mass_tolerance_unit") == "ppm");
    if (!std::isfinite(search_param.precursor_mass_tolerance) || search_param.precursor_mass_tolerance < 0.0 ||
        !std::isfinite(search_param.fragment_mass_tolerance) || search_param.fragment_mass_tolerance < 0.0 ||
        (search_param.fragment_mass_tolerance == 0.0 &&
         (getStringOption_("scoring:method") == "mvh" || getFlag_("scoring:use_mass_accuracy"))))
    {
      OPENMS_LOG_ERROR << "Mass tolerances must be finite and nonnegative; MVH and mass-accuracy weighting require positive fragment tolerance." << endl;
      return ILLEGAL_PARAMETERS;
    }
    search_param.min_length = getIntOption_("oligo:min_size");
    search_param.max_length = getIntOption_("oligo:max_size");
    StringList fixed_mod_names = getStringList_("modifications:fixed");
    search_param.fixed_mods.insert(fixed_mod_names.begin(),
                                   fixed_mod_names.end());

    StringList var_mod_names = getStringList_("modifications:variable");
    search_param.variable_mods.insert(var_mod_names.begin(),
                                      var_mod_names.end());

    resolve_ambiguous_mods_ = getFlag_("modifications:resolve_ambiguities");

    set<ConstRibonucleotidePtr> fixed_modifications =
      getModifications_(search_param.fixed_mods);
    set<ConstRibonucleotidePtr> variable_modifications =
      getModifications_(search_param.variable_mods);
    resolve_ambiguous_mods_ = resolve_ambiguous_mods_ && !ambiguous_mods_.empty();
    RNaseDigestion::CleavageSensitiveModGroups cleavage_sensitive_groups;
    set<ConstRibonucleotidePtr> cleavage_sensitive_modifications;
    set<std::string> cleavage_sensitive_mod_codes;
    set<ConstRibonucleotidePtr> post_digest_variable_modifications =
      variable_modifications;


    // @TODO: add slots for these to "owning run processing metadata"?
    IntList precursor_isotopes = (use_avg_mass ? vector<Int>(1, 0) :
                                  getIntList_("precursor:isotopes"));
    Size max_variable_mods_per_oligo =
      getIntOption_("modifications:variable_max_per_oligo");
    Size report_top_hits = getIntOption_("report:top_hits");

    std::string enzyme_name = getStringOption_("oligo:enzyme");
    const DigestionEnzyme* configured_enzyme =
      RNaseDB::getInstance()->getEnzyme(enzyme_name);
    // the variable modifications that block cuts of @p enzyme (none without an enzyme) are placed during the digestion,
    // the others afterwards
    const auto infer_cleavage_sensitive_mods = [&](const DigestionEnzyme* enzyme)
    {
      cleavage_sensitive_groups = {};
      if (enzyme)
      {
        RNaseDigestion cleavage_sensitive_digestor;
        cleavage_sensitive_digestor.setEnzyme(enzyme);
        cleavage_sensitive_groups =
          cleavage_sensitive_digestor.inferCleavageSensitiveMods(
            variable_modifications);
      }
      cleavage_sensitive_modifications = cleavage_sensitive_groups.combined();
      cleavage_sensitive_mod_codes.clear();
      for (ConstRibonucleotidePtr mod : cleavage_sensitive_modifications)
      {
        cleavage_sensitive_mod_codes.insert(mod->getCode());
      }
      post_digest_variable_modifications = variable_modifications;
      for (ConstRibonucleotidePtr mod : cleavage_sensitive_modifications)
      {
        post_digest_variable_modifications.erase(mod);
      }
    };
    infer_cleavage_sensitive_mods(configured_enzyme);
    if (!cleavage_sensitive_modifications.empty())
    {
      OPENMS_LOG_INFO
        << "Detected " << cleavage_sensitive_modifications.size()
        << " cleavage-sensitive variable modification(s) from enzyme cleavage regexes "
        << "(" << cleavage_sensitive_groups.cuts_before_sensitive.size()
        << " cuts-before-sensitive, "
        << cleavage_sensitive_groups.cuts_after_sensitive.size()
        << " cuts-after-sensitive)."
        << endl;
    }

    StringList potential_adducts =
      getStringList_("precursor:potential_adducts");
    // @TODO: allow different adducts with same mass?
    map<double, std::optional<AdductInfo>> adduct_masses;
    adduct_masses[0.0] = std::nullopt; // always consider "no adduct"
    bool use_adducts = getFlag_("precursor:use_adducts");
    bool include_unknown_charge = getFlag_("precursor:include_unknown_charge");
    bool single_charge_spectra = getFlag_("decharge_ms2");
    if (use_adducts)
    {
      for (const std::string& adduct_name : potential_adducts)
      {
        AdductInfo adduct = parseAdduct_(adduct_name);
        double mass = adduct.getMassShift(use_avg_mass);
        adduct_masses[mass] = adduct;
        OPENMS_LOG_DEBUG << "Added adduct: " << adduct_name << ", mass shift: "
                         << mass << endl;
      }
    }

    // load MS2 map
    MSExperiment spectra;
    FileHandler f;
    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.setOptions(options);
    f.loadExperiment(in_mzml, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW}, log_type_);
    spectra.sortSpectra(true);

    IdentificationData::SourceFile input;
    input.path = test_mode_ ? File::basename(in_mzml) : in_mzml;
    result_run.addSource(input);
    const std::string scoring_method = getStringOption_("scoring:method");
    IdentificationData::ScoreDefinition primary_score;
    primary_score.name = scoring_method == "mvh" ? "MVH score" : "hyperscore";
    primary_score.higher_better = true;
    result_run.setPrimaryScore(result_run.addScore(primary_score));
    IdentificationData::RunSettings settings;
    settings.software = toolName_();
    settings.software_version = test_mode_ ? "test" : version_;
    settings.date = DateTime::now();
    if (test_mode_) settings.date.set("1999-12-31 23:59:59");
    // the raw files behind the mzML input ("spectra_data_raw" in idXML); the mzML itself is the source of the run
    StringList raw_files;
    spectra.getPrimaryMSRunPath(raw_files);
    if (! raw_files.empty()) settings.setMetaValue("spectra_data_raw", raw_files);
    auto& parameters = settings.search;
    parameters.charges = ListUtils::concatenate(search_param.charges, ",");
    parameters.mass_type = use_avg_mass ? SearchParameters::PeakMassType::AVERAGE : SearchParameters::PeakMassType::MONOISOTOPIC;
    parameters.precursor_mass_tolerance = search_param.precursor_mass_tolerance;
    parameters.precursor_mass_tolerance_ppm = search_param.precursor_tolerance_ppm;
    parameters.fragment_mass_tolerance = search_param.fragment_mass_tolerance;
    parameters.fragment_mass_tolerance_ppm = search_param.fragment_tolerance_ppm;
    parameters.fixed_modifications.assign(search_param.fixed_mods.begin(), search_param.fixed_mods.end());
    parameters.variable_modifications.assign(search_param.variable_mods.begin(), search_param.variable_mods.end());
    parameters.setMetaValue("rna_min_length", static_cast<unsigned int>(search_param.min_length));
    parameters.setMetaValue("rna_max_length", static_cast<unsigned int>(search_param.max_length));
    result_run.setSettings(settings);

    // get digested sequences:
    std::string decoy_pattern = getStringOption_("fdr:decoy_pattern");
    if (in_digest.empty()) // new digestion
    {
      search_param.database = in_db;
      search_param.missed_cleavages = getIntOption_("oligo:missed_cleavages");
      search_param.digestion_enzyme = configured_enzyme;
      parameters.missed_cleavages = search_param.missed_cleavages;
      parameters.setMetaValue("rna_enzyme", enzyme_name);
      result_run.setSettings(settings);

      RNaseDigestion digestor;
      digestor.setEnzyme(search_param.digestion_enzyme);
      digestor.setMissedCleavages(search_param.missed_cleavages);
      // set minimum and maximum size of oligo after digestion
      Size min_oligo_length = getIntOption_("oligo:min_size");
      Size max_oligo_length = getIntOption_("oligo:max_size");

      progresslogger.startProgress(0, 1, "loading database from FASTA file...");
      vector<FASTAFile::FASTAEntry> fasta_db;
      FASTAFile().load(in_db, fasta_db);
      progresslogger.endProgress();

      OPENMS_LOG_INFO << "Performing in-silico digestion..." << endl;
      IdentificationData::Database database;
      database.path = in_db;
      IdentificationDataConverter::importSequences(result_run, database, fasta_db, decoy_pattern);
      // Fixed internal modifications affect RNase cleavage. Terminal fixed
      // modifications are applied to each resulting oligo during candidate expansion.
      set<ConstRibonucleotidePtr> internal_fixed;
      for (auto* mod : fixed_modifications)
      {
        if (mod->getTermSpecificity() == Ribonucleotide::ANYWHERE)
        {
          internal_fixed.insert(mod);
        }
      }
      // one candidate per oligo, with its occurrences in the database sequences (as RNaseDigestion::digest for a run)
      map<NASequence, RNaseDigestion::DigestedOligo> candidates;
      if (result_run.getDatabaseSequences())
      {
        for (const auto& entry : *result_run.getDatabaseSequences())
        {
          const NASequence rna = NASequence::fromString(entry.sequence);
          NASequence parent = rna;
          ModifiedNASequenceGenerator::applyFixedModifications(internal_fixed, parent);
          vector<RNaseDigestion::DigestionProduct> products;
          digestor.digestWithCleavageSensitiveMods(parent, cleavage_sensitive_groups,
            max_variable_mods_per_oligo, products, min_oligo_length, max_oligo_length);
          for (const auto& product : products)
          {
            auto [it, inserted] = candidates.try_emplace(product.fragment);
            if (inserted) it->second.sequence = product.fragment;
            // the evidence describes the database sequence, without the modifications of the search
            it->second.addEvidence(entry, rna, product.position);
          }
        }
      }
      digest.reserve(candidates.size());
      for (auto& [sequence, candidate] : candidates)
      {
        digest.push_back(std::move(candidate));
      }

      std::string digest_out = getStringOption_("digest_out");
      if (!digest_out.empty())
      {
        IdentificationData catalog;
        auto& catalog_run = catalog.addRun("digestion", IdentificationData::MoleculeKind::OLIGONUCLEOTIDE);
        auto catalog_settings = settings;
        catalog_settings.setMetaValue("identification:catalog", "true");
        // the catalog comes from the database, not from the spectra
        catalog_settings.removeMetaValue("spectra_data_raw");
        catalog_run.setSettings(catalog_settings);
        // the catalog's databases in the same order, so that the evidence of the candidates refers to them
        for (const auto& database : result_run.getDatabases())
          catalog_run.addDatabase(database);
        catalog_run.setDatabaseSequences(result_run.getDatabaseSequences());
        IdentificationData::SourceFile database_source;
        database_source.path = in_db;
        database_source.identifier = "in-silico digestion";
        auto source = catalog_run.addSource(database_source);
        Size index = 0;
        for (const auto& candidate : digest)
        {
          IdentificationData::Observation observation;
          observation.data_id = "digest=" + std::to_string(++index);
          auto query = catalog_run.addIdentification(source, observation);
          IdentificationData::MatchData value;
          value.encoding = IdentificationData::Encoding::NA_SEQUENCE;
          value.representation = candidate.sequence.toString();
          value.sequence_evidence = candidate.sequence_evidence;
          value.target_decoy = candidate.target_decoy;
          catalog_run.addMatch(query, value);
        }
        FileHandler().storeIdentifications(digest_out, catalog, {FileTypes::IDPARQUET}, log_type_);
      }
    }
    else // load digestion results from a previous run
    {
      OPENMS_LOG_INFO << "Loading pre-digested sequence data..." << endl;
      IdentificationData catalog;
      FileHandler().loadIdentifications(in_digest, catalog, {FileTypes::IDPARQUET}, log_type_);
      if (catalog.getRuns().size() != 1 || catalog.getRuns().front().getMoleculeKind() != IdentificationData::MoleculeKind::OLIGONUCLEOTIDE)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Digest input requires one oligonucleotide catalog run");
      const auto& catalog_run = catalog.getRuns().front();
      // the databases in the same order, so that the evidence of the candidates refers to them
      for (const auto& database : catalog_run.getDatabases())
        result_run.addDatabase(database);
      result_run.setDatabaseSequences(catalog_run.getDatabaseSequences());
      // the digestion settings (missed cleavages, enzyme, oligo lengths) come with the digest; tolerances and
      // modifications are those of this search
      const auto& digestion = catalog_run.getSettings().search;
      parameters.missed_cleavages = digestion.missed_cleavages;
      parameters.enzyme_term_specificity = digestion.enzyme_term_specificity;
      for (const char* key : {"rna_enzyme", "rna_min_length", "rna_max_length"})
        if (digestion.metaValueExists(key)) parameters.setMetaValue(key, digestion.getMetaValue(key));
      // only the enzyme of the digest (if it names one) classifies cleavage-sensitive modifications
      search_param.digestion_enzyme = digestion.metaValueExists("rna_enzyme") ?
        RNaseDB::getInstance()->getEnzyme(digestion.getMetaValue("rna_enzyme").toString()) : nullptr;
      infer_cleavage_sensitive_mods(search_param.digestion_enzyme);
      result_run.setSettings(settings);
      for (const auto& source : catalog_run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            digest.push_back({NASequence::fromString(match.representation), match.sequence_evidence, match.target_decoy});
    }
    Size n_nucleic_acids = result_run.getDatabaseSequences() ? result_run.getDatabaseSequences()->size() : 0;
    if (! decoy_pattern.empty()
        && (! result_run.getDatabaseSequences()
            || std::none_of(result_run.getDatabaseSequences()->begin(), result_run.getDatabaseSequences()->end(), [](const auto& parent) {
             return parent.target_decoy == IdentificationData::TargetDecoy::DECOY;
           })))
    {
      OPENMS_LOG_ERROR << "Error: decoy pattern is set, but no decoy sequences were found" << endl;
      return ILLEGAL_PARAMETERS;
    }

    progresslogger.startProgress(0, 1, "filtering spectra...");
    // @TODO: move this into the loop below (run only when checks pass)
    const std::string deisotoping = getStringOption_("preprocessing:deisotoping");
    bool use_window_mower = getFlag_("preprocessing:filter_window_mower");
    double window_size = getDoubleOption_("preprocessing:window_mower:windowsize");
    int window_peakcount = getIntOption_("preprocessing:window_mower:peakcount");
    std::string window_movetype = getStringOption_("preprocessing:window_mower:movetype");
    bool use_nlargest = getFlag_("preprocessing:filter_nlargest");
    int nlargest_n = getIntOption_("preprocessing:nlargest:n");
    bool remove_precursor = getFlag_("preprocessing:remove_precursor_peak");
    double precursor_mass_tolerance = getDoubleOption_("preprocessing:precursor_mass_tolerance");
    std::string precursor_tolerance_unit = getStringOption_("preprocessing:precursor_mass_tolerance_unit");
    bool precursor_tolerance_ppm = (precursor_tolerance_unit == "ppm");
    int precursor_peak_isotopes = getIntOption_("preprocessing:precursor_peak_isotopes");
    preprocessSpectra_(spectra, search_param.fragment_mass_tolerance,
                       search_param.fragment_tolerance_ppm,
                       single_charge_spectra, negative_mode, min_charge,
                       max_charge, include_unknown_charge,
                       use_window_mower, window_size, window_peakcount, window_movetype,
                       use_nlargest, nlargest_n,
                       remove_precursor, precursor_mass_tolerance,
                       precursor_tolerance_ppm, precursor_peak_isotopes,
                       deisotoping);
    const std::string out_preprocessed = getStringOption_("out_preprocessed");
    if (!out_preprocessed.empty())
    {
      FileHandler().storeExperiment(out_preprocessed, spectra, {FileTypes::MZML}, log_type_);
    }
    progresslogger.endProgress();
    OPENMS_LOG_DEBUG << "preprocessed spectra: " << spectra.getNrSpectra()
                     << endl;

    // build multimap of precursor mass to scan index (and other information):
    multimap<double, PrecursorInfo> precursor_mass_map;
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end();
         ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      const vector<Precursor>& precursors = s_it->getPrecursors();

      // there should be only one precursor and MS2 should contain at least a
      // few peaks to be considered (at least one per nucleotide in the chain):
      if ((precursors.size() != 1) || (s_it->size() < search_param.min_length))
      {
        continue;
      }

      set<Int> possible_charges;
      Int precursor_charge = abs(precursors[0].getCharge());
      if (precursor_charge == 0) // charge information missing
      {
        if (include_unknown_charge)
        {
          possible_charges = search_param.charges; // try all allowed charges
        }
        else
        {
          continue; // skip
        }
      }
      // compare to charge parameters (the charge value in mzML seems to be
      // always positive, so compare by absolute value in negative mode):
      else if ((negative_mode &&
                ((precursor_charge > abs(*search_param.charges.begin())) ||
                 (precursor_charge < abs(*(--search_param.charges.end()))))) ||
               (!negative_mode &&
                ((precursor_charge < *search_param.charges.begin()) ||
                 (precursor_charge > *(--search_param.charges.end())))))
      {
        continue; // charge not in user-supplied range
      }
      else
      {
        possible_charges.insert(precursor_charge); // only one possibility
      }

      double precursor_mz = precursors[0].getMZ();

      for (Int precursor_charge : possible_charges)
      {
        precursor_charge = abs(precursor_charge); // adjust for neg. mode

        // calculate precursor mass (optionally corrected for adducts and peak
        // misassignment) and map it to MS scan index etc.:
        for (const auto& adduct_pair : adduct_masses)
        {
          for (Int isotope_number : precursor_isotopes)
          {
            double precursor_mass =
              calculatePrecursorMass_(precursor_mz, precursor_charge,
                                      isotope_number, adduct_pair.first,
                                      negative_mode);
            PrecursorInfo info(scan_index, precursor_charge, isotope_number,
                               adduct_pair.second);
            precursor_mass_map.insert(make_pair(precursor_mass, info));
          }
        }
      }
    }

    // create spectrum generator:
    NucleicAcidSpectrumGenerator spectrum_generator;
    Param param = spectrum_generator.getParameters();
    vector<std::string> temp = getStringList_("fragment:ions");
    set<std::string> selected_ions(temp.begin(), temp.end());
    if (resolve_ambiguous_mods_ && !selected_ions.contains("a-B"))
    {
      OPENMS_LOG_WARN << "Warning: option 'modifications:resolve_ambiguities' requires a-B ions in parameter 'fragment:ions' - disabling the option." << endl;
      resolve_ambiguous_mods_ = false;
    }
    for (const auto& code : fragment_ion_codes_)
    {
      std::string param_name = "add_" + code + "_ions";
      if (selected_ions.contains(code))
      {
        param.setValue(param_name, "true");
      }
      else
      {
        param.setValue(param_name, "false");
      }
    }
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    param.setValue("add_precursor_peaks", "false");
    spectrum_generator.setParameters(param);

    vector<HitsByScore> annotated_hits(spectra.size());
    MSExperiment exp_ms2_spectra, theo_ms2_spectra; // debug output

    std::string msg = "scoring oligonucleotide models against spectra...";
    progresslogger.startProgress(0, digest.size(), msg);
    Size hit_counter = 0;

    Int base_charge = negative_mode ? -1 : 1;

// shorter oligos take (possibly much) less time to process than longer ones;
// due to the sorting order of "NASequence", they also appear earlier in the
// container - therefore use dynamic scheduling to distribute work evenly:
#pragma omp parallel for schedule(dynamic)
    for (SignedSize index = 0; index < SignedSize(digest.size()); ++index)
    {
      IF_MASTERTHREAD
      {
        progresslogger.setProgress(index);
      }

      const RNaseDigestion::DigestedOligo* oligo_ref = &digest[index];
      vector<NASequence> all_modified_oligos;
      NASequence ns = oligo_ref->sequence;
      // Only a known digest enzyme can classify existing modifications as
      // cleavage-sensitive placements that consume the variable-modification budget.
      Size consumed_cleavage_sensitive_mods = search_param.digestion_enzyme ?
        countCleavageSensitiveMods_(ns, cleavage_sensitive_mod_codes) : 0;
      if (consumed_cleavage_sensitive_mods > max_variable_mods_per_oligo)
      {
        continue;
      }
      Size remaining_variable_mods =
        max_variable_mods_per_oligo - consumed_cleavage_sensitive_mods;
      ModifiedNASequenceGenerator::applyFixedModifications(
        fixed_modifications, ns);
      ModifiedNASequenceGenerator::applyVariableModifications(
        post_digest_variable_modifications, ns, remaining_variable_mods,
        all_modified_oligos, true);
      OPENMS_LOG_DEBUG << "Oligo: " << ns.toString() << endl;

      // group modified oligos by precursor mass - oligos with the same
      // combination of mods (just different placements) will have same mass:
      map<double, vector<const NASequence*>> modified_oligos_by_mass;
      for (const NASequence& seq : all_modified_oligos)
      {
        double mass = (use_avg_mass ? seq.getAverageWeight() :
                       seq.getMonoWeight());
        modified_oligos_by_mass[mass].push_back(&seq);
      }

      for (const auto& pair : modified_oligos_by_mass)
      {
        double candidate_mass = pair.first;

        // determine MS2 precursors that match to the current mass:
        double tol = search_param.precursor_mass_tolerance;
        if (search_param.precursor_tolerance_ppm)
        {
          tol *= candidate_mass * 1e-6;
        }
        multimap<double, PrecursorInfo>::const_iterator low_it =
          precursor_mass_map.lower_bound(candidate_mass - tol), up_it =
          precursor_mass_map.upper_bound(candidate_mass + tol);

        if (low_it == up_it) continue; // no matching precursor in data

        // collect all relevant charge states for theoret. spectrum generation:
        set<Int> precursor_charges;
        for (auto prec_it = low_it; prec_it != up_it; ++prec_it) // OMS_CODING_TEST_EXCLUDE
        {
          precursor_charges.insert(prec_it->second.charge * base_charge);
        }

        for (const NASequence* seq_ptr : pair.second)
        {
          const NASequence& candidate = *seq_ptr;
          OPENMS_LOG_DEBUG << "Candidate: " << candidate.toString() << " ("
                           << float(candidate_mass) << " Da)" << endl;

          // pre-generate spectra:
          map<Int, MSSpectrum> theo_spectra_by_charge;
          spectrum_generator.getMultipleSpectra(theo_spectra_by_charge,
                                                candidate, precursor_charges,
                                                base_charge);

          const auto blacklist_mz = getDoubleList_("scoring:blacklist_mz");
          if (!blacklist_mz.empty())
          {
            for (auto& [charge, spectrum] : theo_spectra_by_charge)
            {
              filterBlacklistedIons_(spectrum, blacklist_mz, getDoubleOption_("scoring:blacklist_tolerance"));
            }
          }

          for (auto prec_it = low_it; prec_it != up_it; ++prec_it) // OMS_CODING_TEST_EXCLUDE
          {
            OPENMS_LOG_DEBUG << "Matching precursor mass: "
                             << float(prec_it->first) << endl;

            Size charge = prec_it->second.charge;
            // look up theoretical spectrum for this charge:
            MSSpectrum& theo_spectrum =
              theo_spectra_by_charge[charge * base_charge];

            Size scan_index = prec_it->second.scan_index;
            const MSSpectrum& exp_spectrum = spectra[scan_index];
            vector<PeptideHit::PeakAnnotation> annotations;
            // Compute score using selected method
            double score;
            if (scoring_method == "mvh")
            {
              Size num_intensity_classes = getIntOption_("scoring:num_intensity_classes");
              double tic_fraction = getDoubleOption_("scoring:tic_fraction");
              bool use_mass_accuracy = getFlag_("scoring:use_mass_accuracy");
              bool single_peak_credit = getStringOption_("scoring:single_peak_credit") == "true";
              score = MetaboliteSpectralMatching::computeMVHScore(
                search_param.fragment_mass_tolerance,
                search_param.fragment_tolerance_ppm, exp_spectrum, theo_spectrum,
                annotations, num_intensity_classes, tic_fraction, 0.0, use_mass_accuracy, single_peak_credit);
            }
            else // hyperscore
            {
              bool use_mass_accuracy = getFlag_("scoring:use_mass_accuracy");
              score = MetaboliteSpectralMatching::computeHyperScore(
                search_param.fragment_mass_tolerance,
                search_param.fragment_tolerance_ppm, exp_spectrum, theo_spectrum,
                annotations, 0.0, use_mass_accuracy);
            }

            if (!exp_ms2_out.empty())
            {
#pragma omp critical (exp_ms2_out)
              exp_ms2_spectra.addSpectrum(exp_spectrum);
            }
            if (!theo_ms2_out.empty())
            {
              theo_spectrum.setName(candidate.toString());
              theo_spectrum.setRT(exp_spectrum.getRT());
#pragma omp critical (theo_ms2_out)
              theo_ms2_spectra.addSpectrum(theo_spectrum);
            }

            if (!std::isfinite(score) || score < 1e-16)
            {
              continue; // no finite hit
            }

            // Check fragment coverage if required
            double min_coverage = getDoubleOption_("report:require_coverage");
            if (min_coverage > 0.0)
            {
              double actual_coverage = calculateFragmentCoverage_(candidate, annotations);
              if (actual_coverage < min_coverage)
              {
                OPENMS_LOG_DEBUG << "Skipping hit due to insufficient fragment coverage: "
                                 << candidate.toString() << " (coverage: " << actual_coverage
                                 << ", required: " << min_coverage << ")" << endl;
                continue; // Insufficient fragment coverage
              }
            }

#pragma omp atomic
            ++hit_counter;

            OPENMS_LOG_DEBUG << "Score: " << score << endl;

#pragma omp critical (annotated_hits_access)
            {
              HitsByScore& scan_hits = annotated_hits[scan_index];
              if (report_top_hits == 0 || scan_hits.size() < report_top_hits || score >= scan_hits.rbegin()->first)
              {
                AnnotatedHit new_hit;
                new_hit.oligo_ref = oligo_ref;
                new_hit.sequence = candidate;
                new_hit.precursor_error_ppm =
                  (prec_it->first - candidate_mass) / candidate_mass * 1.0e6;
                new_hit.annotations = annotations;
                new_hit.precursor_ref = &(prec_it->second);
                retainTopHit_(scan_hits, score, std::move(new_hit), report_top_hits);
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    OPENMS_LOG_INFO << "Undigested nucleic acids: " << n_nucleic_acids << "\nOligonucleotides: " << digest.size()
                    << "\nSearch hits (spectrum matches): " << hit_counter << endl;

    if (!exp_ms2_out.empty())
    {
      FileHandler().storeExperiment(exp_ms2_out, exp_ms2_spectra, {FileTypes::MZML}, log_type_);
    }
    if (!theo_ms2_out.empty())
    {
      FileHandler().storeExperiment(theo_ms2_out, theo_ms2_spectra, {FileTypes::MZML}, log_type_);
    }

    progresslogger.startProgress(0, 1, "post-processing search hits...");
    postProcessHits_(spectra, annotated_hits, id_data, negative_mode);
    progresslogger.endProgress();
    OPENMS_LOG_INFO << "Identified spectra: " << id_data.getRun("NASE").getNumberOfIdentifications() << endl;

    // FDR:
    if (!decoy_pattern.empty())
    {
      OPENMS_LOG_INFO << "Performing FDR calculations..." << endl;
      calculateAndFilterFDR_(id_data, report_top_hits == 1);
    }
    // Report only database sequences with remaining evidence after FDR filtering (the digest catalog
    // keeps all of them), and calculate their coverage from unique intervals.
    auto& coverage_run = id_data.getRun("NASE");
    using SequenceKey = std::pair<UInt32, std::string>; // database and accession
    std::set<SequenceKey> referenced;
    std::map<SequenceKey, std::set<std::pair<UInt64, UInt64>>> covered;
    for (const auto& source : coverage_run.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
          for (const auto& evidence : match.sequence_evidence)
          {
            const SequenceKey key {evidence.database.value, evidence.accession};
            referenced.insert(key);
            if (evidence.start && evidence.end) covered[key].emplace(*evidence.start, *evidence.end);
          }
    if (coverage_run.getDatabaseSequences())
    {
      auto parents = *coverage_run.getDatabaseSequences();
      std::erase_if(parents, [&](const auto& parent) { return ! referenced.contains({parent.database.value, parent.accession}); });
      for (auto& parent : parents)
      {
        const auto length = NASequence::fromString(parent.sequence).size();
        UInt64 count = 0, end = 0;
        bool first = true;
        for (const auto& interval : covered[{parent.database.value, parent.accession}])
        {
          if (interval.second >= length)
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RNA sequence evidence exceeds sequence length",
                                          parent.accession);
          if (first || interval.first > end) count += interval.second - interval.first + 1;
          else if (interval.second > end)
            count += interval.second - end;
          end = std::max(end, interval.second);
          first = false;
        }
        parent.setMetaValue("coverage", length ? static_cast<double>(count) / length : 0.0);
      }
      coverage_run.setDatabaseSequences(std::move(parents));
    }

    // store results
    if (!db_out.empty())
    {
      FileHandler().storeIdentifications(db_out, id_data, {FileTypes::IDPARQUET}, log_type_);
    }

    MzTab results = IdentificationDataConverter::exportMzTab(id_data);
    OPENMS_LOG_DEBUG << "Nucleic acid rows: "
                     << results.getNucleicAcidSectionRows().size()
                     << "\nOligonucleotide rows: "
                     << results.getOligonucleotideSectionRows().size()
                     << "\nOligo-spectrum match rows: "
                     << results.getOSMSectionRows().size() << endl;

    MzTabFile().store(out, results);

    // dummy "peptide" results:
    if (!id_out.empty())
    {
      vector<ProteinIdentification> proteins;
      PeptideIdentificationList peptides;
      IdentificationDataConverter::exportIDs(id_data, proteins, peptides);
      FileHandler().storeIdentifications(id_out, proteins, peptides, {FileTypes::IDXML});
    }

    if (!lfq_out.empty())
    {
      generateLFQInput_(id_data, lfq_out);
    }

    return EXECUTION_OK;
  }

};

#ifndef OPENMS_NASE_CLASS_TEST
int main(int argc, const char** argv)
{
  NucleicAcidSearchEngine tool;
  return tool.main(argc, argv);
}
#endif
