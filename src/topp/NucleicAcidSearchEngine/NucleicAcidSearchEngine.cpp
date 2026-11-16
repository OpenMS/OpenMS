// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Samuel Wein, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/MATH/StatisticFunctions.h> // for "median"

#include <tuple>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

// file types
#include <OpenMS/FORMAT/BedRModFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

// digestion enzymes
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/SYSTEM/File.h>

// ribonucleotides
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include "ModifiedNASequenceGenerator.h"
#include "RNADeisotoper.h"
#include "RNAPrecursorCorrector.h"
#include "RNASpectrumGraph.h"
#include "RNATagCandidateFilter.h"
#include "RNAModificationDiagnostics.h"
#include "RNARouteScorer.h"
#include "RNARouteCandidateGenerator.h"
#include "RNARescoringFeatures.h"
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
#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>


#include <algorithm>
#include <atomic>
#include <iostream>
#include <fstream>
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

Experimental RNA mass-tag filtering is available with @p candidates:method=graph.
Use @p candidates:method=graph_audit to measure tag support while retaining all
enumerated candidates and the existing scoring. See @ref TOPP_NASE_GraphSearch
for supported ion series, bounds, fallback, spectrum-specific route generation,
optional route/scaffold scoring, and the FLASH MS1 / legacy MS2 workflow.

With @p rescoring:enable, the best hit of every spectrum is rescored with the
in-process Percolator, and the q-values are calculated on the Percolator score.
See @ref TOPP_NASE_Rescoring for the features, fallback and limitations.

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
    setValidFormats_("digest", {"oms"});

    registerOutputFile_("out", "<file>", "", "Output file: mzTab");
    setValidFormats_("out", ListUtils::create<std::string>("mzTab"));

    registerOutputFile_("id_out", "<file>", "", "Output file: idXML (for visualization in TOPPView)", false);
    setValidFormats_("id_out", ListUtils::create<std::string>("idXML"));

    registerOutputFile_("db_out", "<file>", "", "Output file: oms (SQLite database)", false);
    setValidFormats_("db_out", ListUtils::create<std::string>("oms"));

    registerOutputFile_("digest_out", "<file>", "", "Output file: sequence database digest. Ignored if 'digest' input is used.", false);
    setValidFormats_("digest_out", {"oms"});

    registerOutputFile_("lfq_out", "<file>", "", "Output file: targets for label-free quantification using FeatureFinderMetaboIdent ('id' input)", false);
    setValidFormats_("lfq_out", vector<std::string>(1, "tsv"));

    registerOutputFile_("bedrmod_out", "<file>", "", "Output file: bedRMod v2 RNA modification track", false);
    setValidFormats_("bedrmod_out", {"bedrmod"});

    registerInputFile_("bedrmod_chebi_mapping", "<file>", "", "Optional CSV mapping file for bedRMod export ('mod'/'name', 'chebi_id' and/or 'modomics_short_name' columns; required for unsupported custom modification names)", false, true);
    setValidFormats_("bedrmod_chebi_mapping", {"csv"});
    registerTOPPSubsection_("bedrmod", "Reference metadata for bedRMod export");
    for (const std::string key : {"organism", "assembly", "annotation_source", "annotation_version"})
    {
      registerStringOption_("bedrmod:" + key, "<value>", "NA", "Explicit reference " + key + " for bedRMod export (NA when unknown)", false);
    }

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

    registerFlag_("precursor:correct_with_flash", "Infer MS2 monoisotopic precursor m/z using precursor-directed FLASH RNA on the complete parent MS1. Confident corrected or monoisotopic scans search only additional isotope offset 0 around the corrected m/z. The configured 'precursor:isotopes' list is used only for unresolved or ambiguous scans. Requires monoisotopic masses and ppm precursor tolerance.", false);
    registerTOPPSubsection_("precursor:flash", "FLASH RNA parent-MS1 precursor correction options");
    registerDoubleOption_("precursor:flash:mass_tolerance", "<ppm>", 0.0, "MS1 isotope detection and selected-ion matching tolerance in ppm (0 = use precursor:mass_tolerance)", false, true);
    setMinFloat_("precursor:flash:mass_tolerance", 0.0);
    registerDoubleOption_("precursor:flash:min_cos", "<value>", 0.85, "Minimum charge-specific RNA averagine isotope cosine for an MS1 precursor envelope", false, true);
    setMinFloat_("precursor:flash:min_cos", 0.5);
    setMaxFloat_("precursor:flash:min_cos", 1.0);
    registerDoubleOption_("precursor:flash:min_snr", "<value>", 0.25, "Minimum charge-specific SNR for an MS1 precursor envelope", false, true);
    setMinFloat_("precursor:flash:min_snr", 0.25);
    registerOutputFile_("out_precursor_corrections", "<file>", "", "Tab-delimited FLASH precursor correction audit, one row per MS2 (requires precursor:correct_with_flash)", false, true);

    registerTOPPSubsection_("fragment", "Fragment (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment mass tolerance", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", ListUtils::create<std::string>("Da,ppm"));

    registerStringList_("fragment:ions", "<choice>", fragment_ion_codes_, "Fragment ions to include in theoretical spectra", false);
    setValidStrings_("fragment:ions", fragment_ion_codes_);

    registerOutputFile_("out_preprocessed", "<file>", "", "Optional MS2 spectra after preprocessing", false, true);
    setValidFormats_("out_preprocessed", {"mzML"});
    registerTOPPSubsection_("preprocessing", "Spectrum preprocessing options");
    registerStringOption_("preprocessing:deisotoping", "<method>", "legacy", "Isotope removal: legacy spacing detector, none, or FLASHDeconv RNA averagine (profile or centroided input; requires ppm fragment tolerance)", false, true);
    setValidStrings_("preprocessing:deisotoping", {"legacy", "none", "flash"});
    registerTOPPSubsection_("preprocessing:legacy", "Legacy isotope-spacing detector options");
    registerIntOption_("preprocessing:legacy:min_isopeaks", "<num>", 3, "Minimum isotope peaks required for a legacy fragment charge assignment",
                       false, true);
    setMinInt_("preprocessing:legacy:min_isopeaks", 2);
    setMaxInt_("preprocessing:legacy:min_isopeaks", 20);
    registerDoubleOption_("preprocessing:legacy:mass_tolerance", "<ppm>", 0.0,
                          "Legacy isotope-spacing tolerance in ppm; zero uses the fragment matching tolerance and unit", false, true);
    setMinFloat_("preprocessing:legacy:mass_tolerance", 0.0);
    registerFlag_("preprocessing:legacy:sum_intensities", "Sum assigned isotope-envelope intensities into each retained monoisotopic peak", true);
    registerIntOption_("preprocessing:legacy:unknown_max_charge", "<num>", 0,
                       "Optional maximum fragment charge for legacy isotope detection when precursor charge is unknown; zero preserves skipping "
                       "this detection. Does not assign a precursor charge.",
                       false, true);
    setMinInt_("preprocessing:legacy:unknown_max_charge", 0);
    setMaxInt_("preprocessing:legacy:unknown_max_charge", 20);
    registerTOPPSubsection_("preprocessing:flash", "FLASHDeconv RNA isotope detection options");
    registerDoubleOption_("preprocessing:flash:min_cos", "<value>", 0.85, "Minimum charge-specific RNA averagine isotope cosine", false, true);
    setMinFloat_("preprocessing:flash:min_cos", 0.5);
    setMaxFloat_("preprocessing:flash:min_cos", 1.0);
    registerDoubleOption_("preprocessing:flash:min_snr", "<value>", 0.25, "Minimum charge-specific isotope envelope SNR", false, true);
    setMinFloat_("preprocessing:flash:min_snr", 0.25);
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
    setValidStrings_("scoring:method", ListUtils::create<std::string>("hyperscore,mvh,route"));
    registerFlag_("scoring:labile_scaffolds",
                  "For graph_routes with MVH/HyperScore, score diagnostic-supported labile scaffolds as alternative full-spectrum models. Retain the "
                  "full modified precursor and export the selected model.",
                  true);
    registerFlag_("scoring:route_base_loss",
                  "Include coherent a-B ladders in route scoring as a separate ion series; graph generation keeps conventional residue edges", true);
    registerDoubleOption_("scoring:route_mass_accuracy_tolerance", "<ppm_or_Da>", 0.0,
                          "Gaussian mass-accuracy weighting scale in fragment-tolerance units; zero uses the fragment matching tolerance. "
                          "This changes intensity weights without changing fragment matches or graph edges.",
                          false, true);
    setMinFloat_("scoring:route_mass_accuracy_tolerance", 0.0);
    registerDoubleOption_("scoring:route_local_radius", "<Da>", 50.0, "Neutral-mass radius of the local background around route nodes", false, true);
    setMinFloat_("scoring:route_local_radius", 0.001);
    registerIntOption_("scoring:route_local_peaks", "<num>", 10, "Strongest local peaks per route node; unioned by observed peak identity", false,
                       true);
    setMinInt_("scoring:route_local_peaks", 1);
    registerIntOption_("scoring:route_min_nodes", "<num>", 2, "Minimum observed nodes in each anchored route", false, true);
    setMinInt_("scoring:route_min_nodes", 1);
    registerTOPPSubsection_("diagnostic", "Spectrum-specific modification filtering from raw MS2 evidence");
    registerFlag_("diagnostic:filter", "Require raw single-nucleotide or explicit diagnostic evidence for variable residue modifications", true);
    registerStringOption_("diagnostic:filter_scope", "<scope>", "all",
                          "all requires every variable modification diagnostic; explicit gates only residues with fixed diagnostic ions, "
                          "allowing intact graph-supported modifications without a terminal single-nucleotide ion",
                          false, true);
    setValidStrings_("diagnostic:filter_scope", {"all", "explicit"});
    registerStringList_(
      "diagnostic:fixed_ions", "<code=mz>", {},
      "Explicit diagnostic m/z values for labile modifications, e.g. t6A=279.08418; may repeat a residue code. Values are polarity-specific.", false,
      true);
    registerDoubleOption_("diagnostic:min_relative_intensity", "<fraction>", 0.0, "Minimum diagnostic intensity relative to the raw MS2 base peak",
                          false, true);
    setMinFloat_("diagnostic:min_relative_intensity", 0.0);
    setMaxFloat_("diagnostic:min_relative_intensity", 1.0);
    registerIntOption_("diagnostic:max_charge", "<num>", 1,
                       "Largest absolute charge for generic single-nucleotide diagnostics; charges above one require a raw isotope partner. "
                       "Explicit diagnostic m/z rules remain unchanged.",
                       false, true);
    setMinInt_("diagnostic:max_charge", 1);
    setMaxInt_("diagnostic:max_charge", 6);
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

    registerTOPPSubsection_("candidates", "Candidate selection options");
    registerStringOption_("candidates:method", "<method>", "enumeration",
                          "enumeration; graph_audit (score all variants); graph (filter enumerated variants); graph_generate "
                          "(constrain modification placement); graph_generate_audit (compare placement recall with exhaustive graph filtering); "
                          "graph_routes (generate spectrum-specific placements from anchored database-compatible routes)",
                          false, true);
    setValidStrings_("candidates:method", {"enumeration", "graph", "graph_audit", "graph_generate", "graph_generate_audit", "graph_routes"});
    registerTOPPSubsection_("routes", "Anchored graph placement (graph_routes) and route-scoring bounds; explicit limits fail the run");
    registerStringOption_(
      "routes:completion", "<mode>", "bounded",
      "Uncovered positions: labile adds diagnostic-supported labile modifications; bounded_labile limits their unsupported placements; "
      "bounded also permits a limited number of "
      "nonlabile modifications; all enumerates all permitted "
      "modifications under explicit resource limits",
      false, true);
    setValidStrings_("routes:completion", {"labile", "bounded_labile", "bounded", "all"});
    registerIntOption_("routes:max_unanchored_modifications", "<num>", 2,
                       "Maximum variable modifications without intact modified graph placements in bounded completion; graph-placed "
                       "modifications still use the full variable modification budget",
                       false, true);
    setMinInt_("routes:max_unanchored_modifications", 0);
    registerStringOption_(
      "routes:precursor_index", "<choice>", "true",
      "Index disjoint route pairs and attainable completion masses with conservative intervals; false retains the exhaustive oracle", false, true);
    setValidStrings_("routes:precursor_index", {"true", "false"});
    registerStringOption_("routes:allow_full_length", "<choice>", "true",
                          "Permit a full-length mass route to constrain terminal placements; precursor nodes remain excluded from scoring", false,
                          true);
    setValidStrings_("routes:allow_full_length", {"true", "false"});
    registerStringOption_("routes:terminal_states", "<choice>", "enzyme",
                          "enzyme permits the terminal groups produced by the configured RNase at each parent position, as in enumeration; all "
                          "also tries every 5' OH/phosphate and 3' OH/phosphate/cyclic-phosphate assignment of each oligonucleotide",
                          false, true);
    setValidStrings_("routes:terminal_states", {"enzyme", "all"});
    registerFlag_("routes:verify_candidates",
                  "Compare every generated candidate multiset with exhaustive route pairing; fail if either search is incomplete or differs", true);
    registerIntOption_("routes:max_length", "<num>", 25, "Maximum route length in residues", false, true);
    setMinInt_("routes:max_length", 2);
    setMaxInt_("routes:max_length", 32);
    for (const auto& [name, value] : vector<pair<std::string, Int>> {{"max_peaks", 200},
                                                                     {"max_states", 1000000},
                                                                     {"max_completion_nodes", 1000000000},
                                                                     {"max_candidates", 5000000},
                                                                     {"max_supports", 20000000},
                                                                     {"min_nodes", 2}})
    {
      const std::string description = name == "max_states"
        ? "Maximum states per route-generation or route-scoring model query; exceeding it fails the search"
        : "Route search setting: " + name;
      registerIntOption_("routes:" + name, "<num>", value, description, false, true);
      setMinInt_("routes:" + name, 1);
    }
    registerTOPPSubsection_("graph", "Bounded RNA mass tags; full-spectrum scoring and native FDR are retained");
    registerFlag_("graph:verify_matches", "Compare every prefilter rejection, prepared graph decision and witness with the reference matcher; fail on any mismatch", true);
    registerIntOption_("graph:min_tag_length", "<num>", 2, "Minimum number of observed ladder steps in a tag; each step joins two matched peaks and spans up to graph:max_step residues", false, true);
    setMinInt_("graph:min_tag_length", 1);
    registerIntOption_("graph:max_tag_length", "<num>", 4, "Maximum residues spanned by a tag (at least graph:min_tag_length)", false, true);
    setMinInt_("graph:max_tag_length", 1);
    setMaxInt_("graph:max_tag_length", 32);
    registerIntOption_("graph:max_step", "<num>", 2, "Maximum residues per graph edge (2 bridges one missing ladder peak)", false, true);
    setMinInt_("graph:max_step", 1);
    setMaxInt_("graph:max_step", 2);
    for (const auto& [name, value] : vector<pair<std::string, Int>> {{"max_peaks", 100},
                                                                     {"max_charge", 3},
                                                                     {"max_edges", 20000},
                                                                     {"max_query_edges", 20000},
                                                                     {"max_query_states", 100000},
                                                                     {"max_placement_nodes", 100000},
                                                                     {"max_total_edges", 50000000}})
    {
      const std::string description = (name == "max_peaks" || name == "max_charge") ? "Graph evidence selection limit: " + name
                                      : name == "max_total_edges"                     ? "Graph cache bound over all spectra; exceeding it aborts the run"
                                                                                      : "Graph exploration/cache bound: " + name
                                                                                          + "; exceeding it uses enumeration fallback";
      registerIntOption_("graph:" + name, "<num>", value, description, false, true);
      setMinInt_("graph:" + name, name == "max_peaks" ? 2 : 1);
    }
    registerOutputFile_("out_graph", "<file>", "",
                        "Tab-delimited per-spectrum graph sizes and enumeration fallback reasons (requires a graph candidate mode)", false, true);
    setValidFormats_("out_graph", {"tsv"});

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
    setMinInt_("modifications:variable_max_per_oligo", 0);
    registerFlag_("modifications:resolve_ambiguities", "Attempt to resolve ambiguous modifications (e.g. 'mA?' for 'mA'/'Am') based on a-B ions.\nThis incurs a performance cost because two modifications have to be considered for each case.\nRequires a-B ions to be enabled in parameter 'fragment:ions'.");

    registerTOPPSubsection_("oligo", "Oligonucleotide (digestion) options (ignored if 'digest' input is used)");
    registerIntOption_("oligo:min_size", "<num>", 5, "Minimum size an oligonucleotide must have after digestion to be considered in the search", false);
    registerIntOption_("oligo:max_size", "<num>", 0, "Maximum size an oligonucleotide must have after digestion to be considered in the search, leave at 0 for no limit", false);

    registerIntOption_("oligo:missed_cleavages", "<num>", 1, "Number of missed cleavages", false, false);
    registerStringOption_("oligo:cleavage_model", "<model>", "modification_aware",
                          "modification_aware preserves modification-dependent RNase cleavage; canonical reproduces digestion in unmodified-sequence "
                          "space (graph_routes only)",
                          false, true);
    setValidStrings_("oligo:cleavage_model", {"modification_aware", "canonical"});

    StringList all_enzymes;
    RNaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("oligo:enzyme", "<choice>", "no cleavage", "The enzyme used for RNA digestion", false);
    setValidStrings_("oligo:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Number of top-scoring hits per spectrum to retain; all exact ties at the cutoff are included (0 for all hits)", false, true);
    setMinInt_("report:top_hits", 0);

    registerFlag_("report:precursor_percentage", "Report the strongest peak in the configured precursor m/z window as a percentage of the raw MS2 base peak intensity. Uses its independent report tolerance; wide windows can include isotope peaks. This does not integrate an isotope envelope or measure preprocessing removal, and is omitted unless there is exactly one precursor.", false);
    registerDoubleOption_("report:precursor_percentage_mass_tolerance", "<tolerance>", 10.0, "Mass tolerance for finding the residual precursor peak in the raw MS2 spectrum", false, true);
    setMinFloat_("report:precursor_percentage_mass_tolerance", 0.0);
    registerStringOption_("report:precursor_percentage_mass_tolerance_unit", "<unit>", "ppm", "Unit for the residual precursor peak mass tolerance", false, true);
    setValidStrings_("report:precursor_percentage_mass_tolerance_unit", {"Da", "ppm"});
    registerDoubleOption_("report:require_coverage", "<fraction>", 0.0, "Minimum fraction of internal oligonucleotide positions that must be covered by fragment ions (0.0 = no requirement, 1.0 = full coverage required). Coverage requires at least one fragment of type a, a-B, b, c, d, w, x, y, or z at each position.", false, true);
    setMinFloat_("report:require_coverage", 0.0);
    setMaxFloat_("report:require_coverage", 1.0);

    registerTOPPSubsection_("fdr", "False Discovery Rate options");
    registerStringOption_("fdr:decoy_pattern", "<string>", "", "std::string used as part of the accession to annotate decoy sequences (e.g. 'DECOY_'). Leave empty to skip the FDR/q-value calculation.", false);
    registerDoubleOption_("fdr:cutoff", "<value>", 1.0, "Cut-off for FDR filtering; search hits with higher q-values will be removed", false);
    setMinFloat_("fdr:cutoff", 0.0);
    setMaxFloat_("fdr:cutoff", 1.0);
    registerFlag_("fdr:remove_decoys", "Do not score hits to decoy sequences and remove them when filtering");

    registerTOPPSubsection_("rescoring", "Semi-supervised rescoring of the search hits with in-process Percolator (requires 'fdr:decoy_pattern')");
    registerFlag_("rescoring:enable",
                  "Rescore the best hit of every spectrum with Percolator on match features (score margins to the other hits, precursor "
                  "error, ion-series coverage, matched intensity and the evidence of the candidate method). The q-values are then "
                  "calculated on the Percolator score; the native score, Percolator score and PEP are reported. The order of the hits of a "
                  "spectrum is unchanged, and hits below the best one get no Percolator score or q-value. If Percolator cannot run (e.g. "
                  "too few decoys) or rejects its learned model, the native score is used with a warning.");
    registerIntOption_("rescoring:top_hits", "<num>", 2,
                       "Hits per spectrum retained for the score-margin features of rescoring (at least 'report:top_hits'; 0 for all)",
                       false, true);
    setMinInt_("rescoring:top_hits", 0);
    registerIntOption_("rescoring:min_decoys", "<num>", 100,
                       "Minimum number of spectra whose best hit is a decoy; with fewer, Percolator's model and PEPs are unreliable and the "
                       "native score is used",
                       false, true);
    setMinInt_("rescoring:min_decoys", 1);
    registerIntOption_("rescoring:seed", "<num>", 1, "Random seed of Percolator's cross-validation folds", false, true);
    registerDoubleOption_("rescoring:train_fdr", "<value>", 0.01,
                          "FDR at which Percolator selects positive training examples and checks the learned model on held-out folds "
                          "(its 'train_fdr' and 'test_fdr'). If a fold has no target at this FDR, the model is rejected and the native "
                          "score is used; a higher value can let small or hard searches train.",
                          false, true);
    setMinFloat_("rescoring:train_fdr", 0.001);
    setMaxFloat_("rescoring:train_fdr", 0.5);
    registerOutputFile_("rescoring:out_features", "<file>", "",
                        "Optional output file: the feature row, label and Percolator result of every rescored hit (for inspection)", false, true);
    setValidFormats_("rescoring:out_features", {"tsv"});
  }

  // relevant information about an MS2 precursor ion
  struct PrecursorInfo
  {
    Size scan_index;
    Int charge;
    Size isotope;
    IdentificationData::AdductOpt adduct;

    PrecursorInfo(Size scan_index, Int charge, Size isotope,
                  const IdentificationData::AdductOpt& adduct = std::nullopt):
      scan_index(scan_index), charge(charge), isotope(isotope), adduct(adduct)
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
    IdentificationData::IdentifiedOligoRef oligo_ref;
    NASequence sequence;
    double precursor_error_ppm = 0.0; // precursor mass error in ppm
    vector<PeptideHit::PeakAnnotation> annotations; // peak/ion annotations
    const PrecursorInfo* precursor_ref = nullptr; // precursor information
    bool graph_evaluated = false;
    RNATagCandidateFilter::Evidence graph_evidence;
    RNARouteScorer::Result route_evidence;
    std::string native_scoring_scaffold;
    std::optional<IdentificationData::IdentifiedOligo> parent_override;
    std::optional<double> percolator_score; // set by rescoring (best hits of a spectrum)
    std::optional<double> percolator_pep;
  };

  typedef multimap<double, AnnotatedHit, greater<double>> HitsByScore;

  /// Score both enumerated and graph-generated models through the existing OpenMS scorers.
  static double scoreNativeSpectrum_(const std::string& method,
                                     double tolerance,
                                     bool tolerance_ppm,
                                     const MSSpectrum& experimental,
                                     const MSSpectrum& theoretical,
                                     vector<PeptideHit::PeakAnnotation>& annotations,
                                     Size intensity_classes,
                                     double tic_fraction,
                                     bool mass_accuracy,
                                     bool single_peak_credit)
  {
    if (method == "mvh")
    {
      return MetaboliteSpectralMatching::computeMVHScore(tolerance, tolerance_ppm, experimental, theoretical, annotations, intensity_classes,
                                                         tic_fraction, 0.0, mass_accuracy, single_peak_credit);
    }
    OPENMS_PRECONDITION(method == "hyperscore", "A native scoring method is required");
    return MetaboliteSpectralMatching::computeHyperScore(tolerance, tolerance_ppm, experimental, theoretical, annotations, 0.0, mass_accuracy);
  }

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

  /// @brief Register a digestion product with its original parent coordinates and neighbors.
  void registerDigestedOligo_(IdentificationData& id_data,
                             IdentificationData::ParentSequenceRef parent_ref,
                             const NASequence& parent,
                             const RNaseDigestion::DigestionProduct& product) const
  {
    IdentificationData::IdentifiedOligo oligo(product.fragment);
    const Size start = product.position.first;
    const Size end = start + product.position.second;
    IdentificationData::ParentMatch match(start, end - 1);
    match.left_neighbor = start > 0 ? parent[start - 1]->getCode() :
      std::string(1, IdentificationData::ParentMatch::LEFT_TERMINUS);
    match.right_neighbor = end < parent.size() ? parent[end]->getCode() :
      std::string(1, IdentificationData::ParentMatch::RIGHT_TERMINUS);
    oligo.parent_matches[parent_ref].insert(match);
    id_data.registerIdentifiedOligo(oligo);
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
    bool make_single_charged = true,
                                           bool annotate_charge = false,
                                           bool sum_intensities = false)
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

    std::vector<double> feature_intensities(feature_number, 0.0);
    if (sum_intensities)
    {
      for (Size i = 0; i < old_spectrum.size(); ++i)
      {
        if (features[i] >= 0) { feature_intensities[features[i]] += old_spectrum[i].getIntensity(); }
      }
    }
    in.clear(false);
    MSSpectrum::IntegerDataArray charge_labels;
    charge_labels.setName("NASE:legacy_fragment_charge");
    const auto append = [&](const Peak1D& peak, Int charge, Size index) {
      in.push_back(peak);
      if (sum_intensities && charge != 0) { in.back().setIntensity(feature_intensities[features[index]]); }
      if (annotate_charge) { charge_labels.push_back(make_single_charged && charge != 0 ? 1 : std::abs(charge)); }
    };
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
        { append(old_spectrum[i], z, i);
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
          append(p, z, i);
        }
      }
      else
      {
        // keep all unassigned peaks
        if (features[i] < 0)
        {
          append(old_spectrum[i], 0, i);
          continue;
        }

        // convert mono-isotopic peak with charge assigned by deisotoping
        if (z != 0)
        {
          if (!make_single_charged)
          { append(old_spectrum[i], z, i);
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
            append(p, z, i);
          }
        }
      }
    }

    if (annotate_charge) { in.getIntegerDataArrays().push_back(std::move(charge_labels)); }

    in.sortByPosition();
  }


  /// Return no value for spectra without one selected precursor or a positive base peak.
  static optional<double> calculatePrecursorPercentage_(const MSSpectrum& spectrum,
                                                        double mass_tolerance,
                                                        bool tolerance_ppm)
  {
    if (spectrum.getPrecursors().size() != 1) return nullopt;

    double precursor_mz = spectrum.getPrecursors()[0].getMZ();
    if (!isfinite(precursor_mz) || precursor_mz <= 0.0) return nullopt;
    double tolerance = tolerance_ppm ? precursor_mz * mass_tolerance * 1e-6 : mass_tolerance;
    double base_peak_intensity = 0.0;
    double residual_intensity = 0.0;
    for (const Peak1D& peak : spectrum)
    {
      double intensity = peak.getIntensity();
      if (!isfinite(intensity)) continue;
      base_peak_intensity = max(base_peak_intensity, intensity);
      if (abs(peak.getMZ() - precursor_mz) <= tolerance)
      {
        residual_intensity = max(residual_intensity, intensity);
      }
    }
    if (base_peak_intensity <= 0.0) return nullopt;
    return 100.0 * residual_intensity / base_peak_intensity;
  }

  double flashMassUpperBound_(double mz, Int charge, bool negative_mode) const
  {
    // Cover the same neutral-mass range and isotope margin in MS1 and MS2.
    const double carrier = negative_mode ? -Constants::PROTON_MASS_U : Constants::PROTON_MASS_U;
    return (mz - carrier) * std::abs(charge) + 100.0;
  }

  /// Replace this run's derived metric; never reuse a stale input annotation.
  static void updatePrecursorPercentage_(MSSpectrum& spectrum, bool enabled,
                                        double mass_tolerance, bool tolerance_ppm)
  {
    spectrum.removeMetaValue("precursor_percentage");
    if (enabled)
    {
      const auto percentage = calculatePrecursorPercentage_(spectrum, mass_tolerance, tolerance_ppm);
      if (percentage) spectrum.setMetaValue("precursor_percentage", *percentage);
    }
  }

  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra, bool negative_mode, Int min_charge, Int max_charge, bool include_unknown_charge, bool use_window_mower, double window_size, int window_peakcount, const std::string& window_movetype, bool use_nlargest, int nlargest_n, bool remove_precursor, double precursor_mass_tolerance, bool precursor_tolerance_ppm, int precursor_peak_isotopes, const std::string& deisotoping, double flash_min_cos, double flash_min_snr,
                          bool annotate_legacy_charge = false,
                          Size legacy_min_isopeaks = 3,
                          double legacy_mass_tolerance = 0.0,
                          bool legacy_sum_intensities = false,
                          Size legacy_unknown_max_charge = 0)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

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

    std::unique_ptr<RNADeisotoper> flash_template;
    if (deisotoping == "flash")
    {
      double max_mass = 0.0;
      for (const auto& spec : exp)
      {
        if (!spec.getPrecursors().empty() && spec.getPrecursors()[0].getCharge() != 0)
        {
          const auto& prec = spec.getPrecursors()[0];
          max_mass = std::max(max_mass, flashMassUpperBound_(prec.getMZ(), prec.getCharge(), negative_mode));
        }
        else if (include_unknown_charge && !spec.empty())
        {
          max_mass = std::max(max_mass, flashMassUpperBound_(std::max_element(spec.begin(), spec.end(), [](const auto& a, const auto& b) { return a.getMZ() < b.getMZ(); })->getMZ(), max_charge, negative_mode));
        }
      }
      if (max_mass == 0.0) max_mass = 1000.0; // No spectrum supplied a mass bound.
      flash_template = std::make_unique<RNADeisotoper>(fragment_mass_tolerance, std::abs(max_charge), !negative_mode, max_mass, flash_min_cos, flash_min_snr);
    }
#pragma omp parallel reduction(+: n_zero_charge, n_inferred_charge)
    {
    std::unique_ptr<RNADeisotoper> flash;
    if (flash_template) flash = std::make_unique<RNADeisotoper>(*flash_template);
#pragma omp for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size();
         ++exp_index)
    {
      MSSpectrum& spec = exp[exp_index];
      spec.removeMetaValue("legacy_unknown_fragment_charge_bound");
      // Keep the post-ThresholdMower base peak as the normalization divisor.
      // FLASH must see raw intensities; native IDs need not be unique or present.
      double normalization_scale = 0.0;
      for (const auto& peak : spec) normalization_scale = std::max(normalization_scale, double(peak.getIntensity()));
      const auto normalize = [&spec, normalization_scale]()
      {
        if (normalization_scale > 0.0)
        {
          for (auto& peak : spec) peak.setIntensity(peak.getIntensity() / normalization_scale);
        }
      };

      // sort by mz
      spec.sortByPosition();

      if (spec.getPrecursors().empty())
      {
        normalize();
        continue; // this shouldn't happen
      }
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
        // Fragment isotope spacings can establish their own charges while the
        // precursor charge remains unknown. The opt-in bound limits these tests.
        const Int fragment_charge_bound = precursor_charge ? precursor_charge : Int(legacy_unknown_max_charge);
        if (! precursor_charge && fragment_charge_bound) { spec.setMetaValue("legacy_unknown_fragment_charge_bound", fragment_charge_bound); }
        deisotopeAndSingleChargeMSSpectrum_(spec, coef, coef * fragment_charge_bound,
                                            legacy_mass_tolerance > 0 ? legacy_mass_tolerance : fragment_mass_tolerance,
                                            legacy_mass_tolerance > 0 || fragment_mass_tolerance_unit_ppm, false, legacy_min_isopeaks, 20,
                                            single_charge_spectra, annotate_legacy_charge, legacy_sum_intensities);
      }
      else if (flash && (spec.getPrecursors()[0].getCharge() != 0 || include_unknown_charge))
      {
        flash->filter(spec, int(exp_index), single_charge_spectra);
      }

      normalize();

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
    map<IdentificationData::IdentifiedOligoRef,
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
      const std::string adduct = precursor.adduct ? (*precursor.adduct)->getName() : std::string();
      return make_tuple(hit.sequence.toString(), hit.oligo_ref->sequence.toString(),
                        precursor.charge, precursor.isotope, adduct,
                        precursor.adduct ? (*precursor.adduct)->getEmpiricalFormula().toString() : std::string(),
                        precursor.adduct ? (*precursor.adduct)->getCharge() : 0,
                        precursor.adduct ? (*precursor.adduct)->getMolMultiplier() : UInt(0),
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
      trimToTopHits_(hits, top_hits);
    }
  }

  /// @brief Keep the best @p top_hits hits and every exact tie at the cutoff (0 keeps all hits).
  static void trimToTopHits_(HitsByScore& hits, Size top_hits)
  {
    if ((top_hits == 0) || (hits.size() <= top_hits)) return;
    auto cutoff = hits.begin();
    std::advance(cutoff, top_hits - 1);
    hits.erase(hits.upper_bound(cutoff->first), hits.end());
  }

  static bool isDecoyHit_(const AnnotatedHit& hit)
  {
    return (hit.parent_override ? *hit.parent_override : *hit.oligo_ref).allParentsAreDecoys();
  }

  /// Score types of Percolator rescoring (registered only if 'rescoring:enable' is set).
  struct RescoringScores
  {
    IdentificationData::ScoreTypeRef score_ref;
    IdentificationData::ScoreTypeRef pep_ref;
  };

  /**
    @brief Rescore the best hit of every spectrum with Percolator.

    The hits must be sorted and ambiguity-resolved. Percolator gets one row per spectrum: its best native hit,
    described relative to the other retained hits (score margins). This is the usual Percolator input (one
    PSM per spectrum with delta scores), so Percolator's PEPs refer to the same spectrum-level competition as
    NASE's FDR calculation. Hits tied with the best native score share its Percolator score and
    PEP; lower-ranked hits get none (and therefore no q-value). The order of the hits of a spectrum is unchanged.

    @return False (no Percolator scores) if Percolator cannot run (e.g. with too few decoys) or its sanity check
            rejects the learned model.
  */
  bool rescoreHits_(const PeakMap& exp, vector<HitsByScore>& annotated_hits)
  {
    Size min_charge = std::numeric_limits<Size>::max(), max_charge = 0, n_rows = 0, n_decoys = 0;
    for (const auto& hits : annotated_hits)
    {
      if (hits.empty()) continue;
      const AnnotatedHit& best = hits.begin()->second;
      const Size charge = static_cast<Size>(std::abs(best.precursor_ref->charge));
      min_charge = std::min(min_charge, charge);
      max_charge = std::max(max_charge, charge);
      ++n_rows;
      if (isDecoyHit_(best)) ++n_decoys;
    }
    const Size min_decoys = getIntOption_("rescoring:min_decoys");
    if ((n_decoys < min_decoys) || (n_decoys == n_rows))
    {
      OPENMS_LOG_WARN << "Warning: Rescoring needs target spectrum matches and at least " << min_decoys
                      << " decoy ones ('rescoring:min_decoys'); found " << n_rows - n_decoys << " target and " << n_decoys
                      << " decoy best hits. Using the native search score." << endl;
      return false;
    }

    // The charge columns cover the charges of all rows (independent of the labels).
    RNARescoringFeatures::Options options;
    options.ion_series = getStringList_("fragment:ions");
    options.min_charge = std::max<Size>(min_charge, 1);
    options.max_charge = std::max(max_charge, options.min_charge);
    const std::string candidate_method = getStringOption_("candidates:method");
    options.route = getStringOption_("scoring:method") == "route";
    options.graph = (candidate_method != "enumeration") && (candidate_method != "graph_routes");
    options.diagnostics = getFlag_("diagnostic:filter");
    options.flash = getFlag_("precursor:correct_with_flash");
    const RNARescoringFeatures builder(options);

    vector<vector<double>> scan_rows(annotated_hits.size());
#pragma omp parallel for
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      const HitsByScore& hits = annotated_hits[scan_index];
      if (hits.empty()) continue;
      const MSSpectrum& spectrum = exp[scan_index];
      RNARescoringFeatures::SpectrumInput spectrum_input;
      for (const auto& peak : spectrum) spectrum_input.total_intensity += peak.getIntensity();
      if (spectrum.metaValueExists("diagnostic_allowed_mods"))
      {
        for (const auto& code : spectrum.getMetaValue("diagnostic_allowed_mods").toStringList())
        {
          spectrum_input.diagnostic_allowed.insert(code);
        }
      }
      if (spectrum.metaValueExists("flash_precursor_cosine"))
      {
        spectrum_input.flash_cosine = double(spectrum.getMetaValue("flash_precursor_cosine"));
        if (spectrum.metaValueExists("flash_precursor_snr")) spectrum_input.flash_snr = double(spectrum.getMetaValue("flash_precursor_snr"));
      }
      vector<RNARescoringFeatures::HitInput> inputs;
      inputs.reserve(hits.size());
      for (const auto& [score, hit] : hits)
      {
        RNARescoringFeatures::HitInput input;
        input.native_score = score;
        input.sequence = hit.sequence;
        input.precursor_error_ppm = hit.precursor_error_ppm;
        input.isotope_offset = hit.precursor_ref->isotope;
        input.charge = hit.precursor_ref->charge;
        input.adduct = hit.precursor_ref->adduct.has_value();
        input.annotations = &hit.annotations;
        input.route_score = hit.route_evidence.score;
        input.route_nodes = hit.route_evidence.nodes;
        input.route_matched_intensity = hit.route_evidence.matched_intensity;
        input.route_local_intensity = hit.route_evidence.local_intensity;
        input.route_prefix_length = hit.route_evidence.prefix_length;
        input.route_suffix_length = hit.route_evidence.suffix_length;
        input.graph_tag_length = (hit.graph_evaluated && hit.graph_evidence.supported) ? hit.graph_evidence.length : 0;
        inputs.push_back(std::move(input));
      }
      // All retained hits give the margins; only the best one is rescored.
      scan_rows[scan_index] = std::move(builder.computeSpectrum(inputs, spectrum_input).front());
    }

    RescoreInput input;
    input.feature_names = builder.featureNames();
    input.features.reserve(n_rows);
    vector<Size> row_scans;
    for (Size scan_index = 0; scan_index < annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) continue;
      input.features.push_back(std::move(scan_rows[scan_index]));
      input.is_decoy.push_back(isDecoyHit_(annotated_hits[scan_index].begin()->second));
      input.scan_numbers.push_back(static_cast<int>(scan_index));
      row_scans.push_back(scan_index);
    }

    Percolator percolator;
    Param param = percolator.getParameters();
    param.setValue("seed", getIntOption_("rescoring:seed"));
    param.setValue("train_fdr", getDoubleOption_("rescoring:train_fdr"));
    param.setValue("test_fdr", getDoubleOption_("rescoring:train_fdr"));
    param.setValue("num_threads", std::max(getIntOption_("threads"), 1));
    percolator.setParameters(param);
    RescoreOutput output;
    try
    {
      output = percolator.rescore(input);
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "Warning: Percolator rescoring failed (" << n_rows << " spectra, " << n_decoys << " decoy best hits): " << e.what()
                      << "\nUsing the native search score." << endl;
      writeRescoringFeatures_(exp, annotated_hits, input, row_scans, nullptr);
      return false;
    }
    if (percolator.usedDefaultDirection())
    {
      // Percolator then scores with a single feature, which can rank far worse than the native score (e.g. when a
      // high-scoring decoy leaves no target at q < 1% on a held-out fold).
      OPENMS_LOG_WARN << "Warning: Percolator's sanity check rejected the learned model (" << n_rows << " spectra, " << n_decoys
                      << " decoy best hits). Using the native search score." << endl;
      writeRescoringFeatures_(exp, annotated_hits, input, row_scans, nullptr);
      return false;
    }
    writeRescoringFeatures_(exp, annotated_hits, input, row_scans, &output);

    for (Size row = 0; row < n_rows; ++row)
    {
      HitsByScore& hits = annotated_hits[row_scans[row]];
      const double best_score = hits.begin()->first;
      for (auto it = hits.begin(); (it != hits.end()) && (it->first == best_score); ++it)
      {
        it->second.percolator_score = output.scores[row];
        it->second.percolator_pep = output.peps[row];
      }
    }
    OPENMS_LOG_INFO << "Rescored the best hits of " << n_rows << " spectra (" << n_decoys << " decoys) with Percolator." << endl;
    return true;
  }

  /// Write the rescoring features of the best hits (and the Percolator results, if any) to 'rescoring:out_features'.
  void writeRescoringFeatures_(const PeakMap& exp, const vector<HitsByScore>& annotated_hits, const RescoreInput& input,
                               const vector<Size>& row_scans, const RescoreOutput* output)
  {
    const std::string path = getStringOption_("rescoring:out_features");
    if (path.empty()) return;
    SVOutStream out(path, "\t", "_", QuotingMethod::NONE);
    out << "native_id" << "scan_index" << "sequence" << "charge" << "decoy";
    for (const auto& name : input.feature_names) out << name;
    out << "percolator_score" << "percolator_q_value" << "percolator_pep" << nl;
    for (Size row = 0; row < row_scans.size(); ++row)
    {
      const Size scan_index = row_scans[row];
      const AnnotatedHit& hit = annotated_hits[scan_index].begin()->second;
      out << exp[scan_index].getNativeID() << scan_index << hit.sequence.toString() << hit.precursor_ref->charge
          << (input.is_decoy[row] ? 1 : 0);
      for (double value : input.features[row]) out << value;
      if (output) out << output->scores[row] << output->q_values[row] << output->peps[row];
      else out << "" << "" << "";
      out << nl;
    }
  }

  /**
    @brief Turn the search hits into observation matches.

    With @p rescoring, the best hits are also rescored (rescoreHits_()). Each spectrum keeps its best
    @p report_top_hits hits by native score (with exact ties).

    @return Whether the hits were rescored
  */
  bool postProcessHits_(const PeakMap& exp,
                        vector<HitsByScore>& annotated_hits,
                        IdentificationData& id_data,
                        bool negative_mode,
                        IdentificationData::InputFileRef file_ref,
                        IdentificationData::ScoreTypeRef score_ref,
                        Size report_top_hits,
                        const std::optional<RescoringScores>& rescoring)
  {
    // Use the current run references; cached entries may sort before them.
    const std::string graph_method = getStringOption_("candidates:method");

#pragma omp parallel for
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      HitsByScore& hits = annotated_hits[scan_index];
      if (hits.empty()) continue;
      sortHitsByValue_(hits);
      if (resolve_ambiguous_mods_ && (hits.size() > 1))
      {
        resolveAmbiguousMods_(hits);
        // Merged ambiguity codes can change the value order of the surviving hits.
        sortHitsByValue_(hits);
      }
    }

    const bool rescored = rescoring && rescoreHits_(exp, annotated_hits);

// @TODO: change OpenMP schedule from default ("static") to "dynamic"/"guided"?
#pragma omp parallel for
    for (SignedSize scan_index = 0;
         scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) continue;

      const MSSpectrum& spectrum = exp[scan_index];
      IdentificationData::Observation obs(spectrum.getNativeID(), file_ref,
                                          spectrum.getRT(),
                                          spectrum.getPrecursors()[0].getMZ());
      obs.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
      obs.setMetaValue("precursor_intensity",
                         spectrum.getPrecursors()[0].getIntensity());
      if (spectrum.metaValueExists("precursor_percentage"))
      {
        obs.setMetaValue("precursor_percentage", spectrum.getMetaValue("precursor_percentage"));
      }
      IdentificationData::ObservationRef obs_ref;
#pragma omp critical (id_data_access)
      obs_ref = id_data.registerObservation(obs);

      // The search may retain more hits than are reported (search_top_hits).
      trimToTopHits_(annotated_hits[scan_index], report_top_hits);

      // create full oligo hit structure from annotated hits
      for (const auto& pair : annotated_hits[scan_index])
      {
        double score = pair.first;
        const AnnotatedHit& hit = pair.second;
        OPENMS_LOG_DEBUG << "Hit sequence: " << hit.sequence.toString() << endl;

        // transfer parent matches from unmodified oligo:
        IdentificationData::IdentifiedOligo oligo = hit.parent_override ? *hit.parent_override : *hit.oligo_ref;
        oligo.sequence = hit.sequence;
        IdentificationData::IdentifiedOligoRef oligo_ref;
#pragma omp critical (id_data_access)
        oligo_ref = id_data.registerIdentifiedOligo(oligo);

        Int charge = hit.precursor_ref->charge;
        if ((charge > 0) && negative_mode) charge = -charge;
        IdentificationData::ObservationMatch match(oligo_ref, obs_ref, charge);
        match.addScore(score_ref, score, id_data.getCurrentProcessingStep());
        if (hit.percolator_score) match.addScore(rescoring->score_ref, *hit.percolator_score, id_data.getCurrentProcessingStep());
        if (hit.percolator_pep) match.addScore(rescoring->pep_ref, *hit.percolator_pep, id_data.getCurrentProcessingStep());
        match.peak_annotations[id_data.getCurrentProcessingStep()] =
          hit.annotations;
        // @TODO: add a field for this to "IdentificationData::ObservationMatch"?
        match.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
                           hit.precursor_error_ppm);
        match.setMetaValue("isotope_offset", hit.precursor_ref->isotope);
        if (graph_method == "graph_routes") { match.setMetaValue("candidate_method", graph_method); }
        if (getStringOption_("scoring:method") == "route")
        {
          match.setMetaValue("route_score", hit.route_evidence.score);
          match.setMetaValue("route_unique_nodes", static_cast<unsigned int>(hit.route_evidence.nodes));
          match.setMetaValue("route_matched_intensity", hit.route_evidence.matched_intensity);
          match.setMetaValue("route_local_intensity", hit.route_evidence.local_intensity);
          match.setMetaValue("route_prefix_length", static_cast<unsigned int>(hit.route_evidence.prefix_length));
          match.setMetaValue("route_prefix_scaffold", hit.route_evidence.prefix_length && hit.route_evidence.prefix_scaffold.empty()
                                                        ? hit.sequence.toString()
                                                        : hit.route_evidence.prefix_scaffold);
          match.setMetaValue("route_suffix_scaffold", hit.route_evidence.suffix_length && hit.route_evidence.suffix_scaffold.empty()
                                                        ? hit.sequence.toString()
                                                        : hit.route_evidence.suffix_scaffold);
          match.setMetaValue("route_suffix_length", static_cast<unsigned int>(hit.route_evidence.suffix_length));
        }
        if (! hit.native_scoring_scaffold.empty()) { match.setMetaValue("native_scoring_scaffold", hit.native_scoring_scaffold); }
        if (spectrum.metaValueExists("diagnostic_allowed_mods"))
        {
          match.setMetaValue("diagnostic_allowed_mods", spectrum.getMetaValue("diagnostic_allowed_mods"));
          match.setMetaValue("diagnostic_labile_mods", spectrum.getMetaValue("diagnostic_labile_mods"));
          match.setMetaValue("diagnostic_matches", spectrum.getMetaValue("diagnostic_matches"));
          match.setMetaValue("diagnostic_match_charges", spectrum.getMetaValue("diagnostic_match_charges"));
          match.setMetaValue("diagnostic_isotope_matches", spectrum.getMetaValue("diagnostic_isotope_matches"));
        }
        if (spectrum.metaValueExists("legacy_unknown_fragment_charge_bound"))
        {
          match.setMetaValue("legacy_unknown_fragment_charge_bound", spectrum.getMetaValue("legacy_unknown_fragment_charge_bound"));
        }
        if (hit.graph_evaluated)
        {
          match.setMetaValue("graph_candidate_method", graph_method);
          match.setMetaValue("graph_candidate_status", hit.graph_evidence.supported ? "supported"
                                                                                    : (hit.graph_evidence.fallback_reason.empty()
                                                                                         ? "rejected_in_audit"
                                                                                         : "fallback:" + hit.graph_evidence.fallback_reason));
          match.setMetaValue("graph_tag_length", static_cast<unsigned int>(hit.graph_evidence.length));
          match.setMetaValue("graph_ion_series", hit.graph_evidence.ion_series);
          match.setMetaValue("graph_fragment_charge", static_cast<unsigned int>(hit.graph_evidence.charge));
          match.setMetaValue("graph_start_ion", static_cast<unsigned int>(hit.graph_evidence.start_ion));
          if (hit.graph_evidence.supported)
          {
            const bool prefix = std::string("abcd").find(hit.graph_evidence.ion_series) != std::string::npos;
            const Size start = prefix ? hit.graph_evidence.start_ion : hit.sequence.size() - hit.graph_evidence.start_ion - hit.graph_evidence.length;
            match.setMetaValue("graph_tag_sequence", hit.sequence.getSubsequence(start, hit.graph_evidence.length).toString());
          }
        }
        match.adduct_opt = hit.precursor_ref->adduct;
        for (const auto& key : {"flash_precursor_status", "flash_precursor_original_mz",
             "flash_precursor_original_charge", "flash_precursor_parent", "flash_precursor_parent_source",
             "flash_precursor_parent_peaks", "flash_precursor_parent_queries",
             "flash_precursor_mz", "flash_precursor_mono_mass", "flash_precursor_isotope",
             "flash_precursor_cosine", "flash_precursor_snr", "flash_precursor_mono_observed"})
        {
          if (exp[scan_index].metaValueExists(key)) match.setMetaValue(key, exp[scan_index].getMetaValue(key));
        }
#pragma omp critical (id_data_access)
        id_data.registerObservationMatch(match);
      }
    }
    id_data.cleanup();
    return rescored;
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


  /// q-values from target-decoy competition on @p score_ref (the score that ranks the hits)
  void calculateAndFilterFDR_(IdentificationData& id_data, bool only_top_hits, IdentificationData::ScoreTypeRef score_ref)
  {
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getDefaults();
    fdr_params.setValue("use_all_hits", only_top_hits ? "false" : "true");
    bool remove_decoys = getFlag_("fdr:remove_decoys");
    fdr_params.setValue("add_decoy_peptides", remove_decoys ? "false" : "true");
    fdr.setParameters(fdr_params);
    IdentificationData::ScoreTypeRef fdr_ref =
    fdr.applyToObservationMatches(id_data, score_ref);
    double fdr_cutoff = getDoubleOption_("fdr:cutoff");
    if (remove_decoys) // remove references to decoys from shared oligos
    {
      IDFilter::removeDecoys(id_data);
    }
    if (fdr_cutoff < 1.0)
    {
      IDFilter::filterObservationMatchesByScore(id_data, fdr_ref, fdr_cutoff);
      OPENMS_LOG_INFO << "Search hits after FDR filtering: "
                      << id_data.getObservationMatches().size()
                      << "\nIdentified spectra after FDR filtering: "
                      << id_data.getObservations().size() << endl;
    }
  }

  /// @brief Emit null columns for unavailable RNA scores.
  static void completeRNAScoreColumns_(MzTab& table)
  {
    const auto& metadata = table.getMetaData();
    const auto complete = [&](auto& rows, const auto& scores) {
      for (auto& row : rows)
      {
        for (const auto& [score, definition] : scores)
        {
          row.best_search_engine_score.try_emplace(score);
          for (const auto& [run, info] : metadata.ms_run)
          {
            row.search_engine_score_ms_run[score].try_emplace(run);
          }
        }
      }
    };
    // Some core versions serialize only present score-map entries. Mixed target/
    // decoy and tied results can have unavailable parent/alternative q-values.
    auto nucleic_acids = table.getNucleicAcidSectionRows();
    complete(nucleic_acids, metadata.nucleic_acid_search_engine_score);
    table.setNucleicAcidSectionRows(nucleic_acids);
    auto oligos = table.getOligonucleotideSectionRows();
    complete(oligos, metadata.oligonucleotide_search_engine_score);
    table.setOligonucleotideSectionRows(oligos);
    auto matches = table.getOSMSectionRows();
    for (auto& row : matches)
    {
      for (const auto& [score, definition] : metadata.osm_search_engine_score)
      {
        row.search_engine_score.try_emplace(score);
      }
    }
    table.setOSMSectionRows(matches);
  }


  void generateLFQInput_(IdentificationData& id_data, const std::string& out_file)
  {
    using AdductedOligo = pair<NASequence, IdentificationData::AdductOpt>;
    using PrecursorPair = pair<double, double>; // precursor intensity, RT
    // mapping: charge -> list of precursors
    using PrecursorsByCharge = map<Int, vector<PrecursorPair>>;
    map<AdductedOligo, PrecursorsByCharge> rt_info;
    for (const IdentificationData::ObservationMatch& match :
           id_data.getObservationMatches())
    {
      const NASequence& seq =
        match.identified_molecule_var.getIdentifiedOligoRef()->sequence;
      auto key = make_pair(seq, match.adduct_opt);
      double rt = match.observation_ref->rt;
      double prec_int =
        match.observation_ref->getMetaValue("precursor_intensity");
      rt_info[key][match.charge].push_back(make_pair(prec_int, rt));
    }

    SVOutStream tsv(out_file);
    tsv.modifyStrings(false);
    tsv << "CompoundName" << "SumFormula" << "Mass" << "Charge"
        << "RetentionTime" << "RetentionTimeRange" << "IsoDistribution" << endl;
    for (const auto& entry : rt_info)
    {
      std::string name = entry.first.first.toString();
      EmpiricalFormula ef = entry.first.first.getFormula();
      const IdentificationData::AdductOpt& adduct = entry.first.second;
      if (adduct)
      {
        name += "+[" + (*adduct)->getName() + "]";
        ef += (*adduct)->getEmpiricalFormula();
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

    IdentificationData id_data; // container for results

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
    const std::string bedrmod_out = getStringOption_("bedrmod_out");
    const std::string bedrmod_chebi_mapping = getStringOption_("bedrmod_chebi_mapping");
    std::string lfq_out = getStringOption_("lfq_out");
    std::string theo_ms2_out = getStringOption_("theo_ms2_out");
    std::string exp_ms2_out = getStringOption_("exp_ms2_out");
    bool use_avg_mass = getFlag_("precursor:use_avg_mass");
    const std::string candidate_method = getStringOption_("candidates:method");
    const bool route_generation = candidate_method == "graph_routes";
    const bool canonical_cleavage = getStringOption_("oligo:cleavage_model") == "canonical";
    const bool graph_enabled = candidate_method != "enumeration" && ! route_generation;
    const bool graph_generate = candidate_method == "graph_generate" || candidate_method == "graph_generate_audit";
    const bool placement_audit = candidate_method == "graph_generate_audit";
    const bool verify_graph_matches = getFlag_("graph:verify_matches");
    const bool graph_filter = candidate_method == "graph" || graph_generate;
    const Size max_placement_nodes = getIntOption_("graph:max_placement_nodes");
    const std::string graph_out = getStringOption_("out_graph");
    if (! graph_enabled && ((! route_generation && ! graph_out.empty()) || verify_graph_matches))
    {
      OPENMS_LOG_ERROR << "'out_graph' and 'graph:verify_matches' require a graph mode in 'candidates:method'.\n";
      return ILLEGAL_PARAMETERS;
    }
    if (graph_enabled && getIntOption_("graph:min_tag_length") > getIntOption_("graph:max_tag_length"))
    {
      OPENMS_LOG_ERROR << "'graph:min_tag_length' must not exceed 'graph:max_tag_length'.\n";
      return ILLEGAL_PARAMETERS;
    }
    if ((canonical_cleavage && ! route_generation) || (route_generation && use_avg_mass))
    {
      OPENMS_LOG_ERROR << "Canonical cleavage requires graph_routes; graph_routes requires monoisotopic masses.\n";
      return ILLEGAL_PARAMETERS;
    }
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

    IdentificationData::DBSearchParam search_param;
    for (Int charge = min_charge; abs(charge) <= abs(max_charge);
         charge += charge_step)
    {
      search_param.charges.insert(charge);
    }
    search_param.molecule_type = IdentificationData::MoleculeType::RNA;
    search_param.mass_type = (use_avg_mass ?
                              IdentificationData::MassType::AVERAGE :
                              IdentificationData::MassType::MONOISOTOPIC);
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
    if (graph_enabled && search_param.fragment_mass_tolerance == 0.0)
    {
      OPENMS_LOG_ERROR << "Graph candidate modes require a positive 'fragment:mass_tolerance'." << endl;
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


    // @TODO: add slots for these to "IdentificationData::DBSearchParam"?
    IntList precursor_isotopes = (use_avg_mass ? vector<Int>(1, 0) :
                                  getIntList_("precursor:isotopes"));
    Size max_variable_mods_per_oligo =
      getIntOption_("modifications:variable_max_per_oligo");
    Size report_top_hits = getIntOption_("report:top_hits");
    // Hits retained per spectrum during the search; postProcessHits_ trims them to report_top_hits.
    // Rescoring retains runner-up hits for its score-margin features (0 = all hits).
    const bool rescoring_enabled = getFlag_("rescoring:enable");
    Size search_top_hits = report_top_hits;
    if (rescoring_enabled && (search_top_hits != 0))
    {
      const Size rescoring_top_hits = getIntOption_("rescoring:top_hits");
      search_top_hits = (rescoring_top_hits == 0) ? 0 : std::max(search_top_hits, rescoring_top_hits);
    }

    std::string enzyme_name = getStringOption_("oligo:enzyme");
    const DigestionEnzyme* configured_enzyme =
      RNaseDB::getInstance()->getEnzyme(enzyme_name);
    RNaseDigestion cleavage_sensitive_digestor;
    cleavage_sensitive_digestor.setEnzyme(configured_enzyme);
    cleavage_sensitive_groups =
      cleavage_sensitive_digestor.inferCleavageSensitiveMods(
        variable_modifications);
    cleavage_sensitive_modifications = cleavage_sensitive_groups.combined();
    for (ConstRibonucleotidePtr mod : cleavage_sensitive_modifications)
    {
      cleavage_sensitive_mod_codes.insert(mod->getCode());
    }
    post_digest_variable_modifications = variable_modifications;
    for (ConstRibonucleotidePtr mod : cleavage_sensitive_modifications)
    {
      post_digest_variable_modifications.erase(mod);
    }
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
    map<double, IdentificationData::AdductOpt> adduct_masses;
    adduct_masses[0.0] = std::nullopt; // always consider "no adduct"
    bool use_adducts = getFlag_("precursor:use_adducts");
    bool include_unknown_charge = getFlag_("precursor:include_unknown_charge");
    bool single_charge_spectra = getFlag_("decharge_ms2");
    const std::string deisotoping = getStringOption_("preprocessing:deisotoping");
    const double flash_min_cos = getDoubleOption_("preprocessing:flash:min_cos");
    const double flash_min_snr = getDoubleOption_("preprocessing:flash:min_snr");
    if (deisotoping == "none" && single_charge_spectra)
    {
      OPENMS_LOG_ERROR << "Decharging requires isotope charge assignment; select legacy or FLASH deisotoping." << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (deisotoping == "flash" && !search_param.fragment_tolerance_ppm)
    {
      OPENMS_LOG_ERROR << "FLASH deisotoping requires ppm fragment tolerance." << endl;
      return ILLEGAL_PARAMETERS;
    }

    if (use_adducts)
    {
      for (const std::string& adduct_name : potential_adducts)
      {
        AdductInfo adduct = parseAdduct_(adduct_name);
        double mass = adduct.getMassShift(use_avg_mass);
        IdentificationData::AdductRef ref = id_data.registerAdduct(adduct);
        adduct_masses[mass] = ref;
        OPENMS_LOG_DEBUG << "Added adduct: " << adduct_name << ", mass shift: "
                         << mass << endl;
      }
    }

    const bool correct_precursor = getFlag_("precursor:correct_with_flash");
    const std::string precursor_report = getStringOption_("out_precursor_corrections");
    if (correct_precursor && use_avg_mass)
    {
      OPENMS_LOG_ERROR << "FLASH precursor correction requires monoisotopic masses; disable use_avg_mass." << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (correct_precursor && !search_param.precursor_tolerance_ppm)
    {
      OPENMS_LOG_ERROR << "FLASH precursor correction requires ppm precursor tolerance." << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (!correct_precursor && !precursor_report.empty())
    {
      OPENMS_LOG_ERROR << "The precursor correction audit requires precursor:correct_with_flash." << endl;
      return ILLEGAL_PARAMETERS;
    }

    // load MS2 map
    MSExperiment spectra;
    FileHandler f;
    PeakFileOptions options;
    options.clearMSLevels();
    if (correct_precursor) options.addMSLevel(1);
    options.addMSLevel(2);
    f.setOptions(options);
    f.loadExperiment(in_mzml, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW}, log_type_);
    const bool report_precursor_percentage = getFlag_("report:precursor_percentage");
    const double report_precursor_tolerance = getDoubleOption_("report:precursor_percentage_mass_tolerance");
    const bool report_precursor_tolerance_ppm = getStringOption_("report:precursor_percentage_mass_tolerance_unit") == "ppm";
    for (MSSpectrum& spectrum : spectra)
    {
      if (!correct_precursor)
      {
        // FLASH annotations belong to the run that produced them. Disabling
        // correction must use the configured isotope list and fresh provenance.
        std::vector<std::string> keys;
        spectrum.getKeys(keys);
        for (const auto& key : keys)
        {
          if (key.starts_with("flash_precursor_")) spectrum.removeMetaValue(key);
        }
      }
      updatePrecursorPercentage_(spectrum, report_precursor_percentage,
        report_precursor_tolerance, report_precursor_tolerance_ppm);
    }
    if (correct_precursor)
    {
      double tolerance = getDoubleOption_("precursor:flash:mass_tolerance");
      if (tolerance == 0) tolerance = search_param.precursor_mass_tolerance;
      if (tolerance <= 0)
      {
        OPENMS_LOG_ERROR << "FLASH precursor correction requires a positive MS1 ppm tolerance." << endl;
        return ILLEGAL_PARAMETERS;
      }
      const int min_abs = std::min(std::abs(min_charge), std::abs(max_charge));
      const int max_abs = std::max(std::abs(min_charge), std::abs(max_charge));
      double max_mass = 0.0;
      for (const auto& spec : spectra)
      {
        if (spec.getMSLevel() != 2) continue;
        for (const auto& p : spec.getPrecursors())
        {
          const int z = p.getCharge() == 0 ? max_abs : std::abs(p.getCharge());
          max_mass = std::max(max_mass, flashMassUpperBound_(p.getMZ(), z, negative_mode));
        }
      }
      if (max_mass == 0.0) max_mass = 1000.0; // No precursor supplied a mass bound.
      RNAPrecursorCorrector corrector(tolerance, min_abs, max_abs, !negative_mode, max_mass,
        getDoubleOption_("precursor:flash:min_cos"), getDoubleOption_("precursor:flash:min_snr"));
      const auto corrections = corrector.correct(spectra);
      std::map<std::string, Size> counts;
      for (const auto& row : corrections) ++counts[row.status];
      OPENMS_LOG_INFO << "FLASH RNA precursor correction: " << corrections.size() << " MS2 scans checked";
      for (const auto& [status, count] : counts) OPENMS_LOG_INFO << "; " << status << "=" << count;
      OPENMS_LOG_INFO << endl;
      if (!precursor_report.empty()) RNAPrecursorCorrector::storeReport(precursor_report, corrections);
      // MS1 scans have served their purpose; preserve the existing MS2-only
      // search indexing and the experiment's file metadata.
      std::vector<MSSpectrum> ms2;
      for (auto& spec : spectra)
      {
        if (spec.getMSLevel() == 2) ms2.push_back(std::move(spec));
      }
      spectra.setSpectra(std::move(ms2));
    }
    spectra.sortSpectra(true);

    // input file meta data:
    std::string input_name = test_mode_ ? File::basename(in_mzml) : in_mzml;
    IdentificationData::InputFile input(input_name);
    vector<std::string> primary_files;
    spectra.getPrimaryMSRunPath(primary_files);
    input.primary_files.insert(primary_files.begin(), primary_files.end());
    IdentificationData::InputFileRef file_ref =
      id_data.registerInputFile(input);
    // processing software meta data:
    const std::string scoring_method = getStringOption_("scoring:method");
    IdentificationData::ScoreType score(scoring_method == "route" ? "RNA route score" : (scoring_method == "mvh" ? "MVH score" : "hyperscore"), true);
    IdentificationData::ScoreTypeRef search_score_ref =
      id_data.registerScoreType(score);
    CVTerm qvalue("MS:1002354", "PSM-level q-value", "MS");
    score = IdentificationData::ScoreType(qvalue, false);
    IdentificationData::ScoreTypeRef qvalue_ref =
      id_data.registerScoreType(score);
    IdentificationData::ProcessingSoftware software(toolName_(), version_);
    // in test mode just overwrite with a generic version:
    if (test_mode_) software.setVersion("test");
    // @TODO: which should be the "primary" (first) score?
    software.assigned_scores.push_back(search_score_ref);
    software.assigned_scores.push_back(qvalue_ref);
    std::optional<RescoringScores> rescoring_scores;
    if (rescoring_enabled)
    {
      RescoringScores refs;
      refs.score_ref = id_data.registerScoreType(IdentificationData::ScoreType(CVTerm("MS:1001492", "percolator:score", "MS"), true));
      refs.pep_ref = id_data.registerScoreType(IdentificationData::ScoreType(CVTerm("MS:1001493", "percolator:PEP", "MS"), false));
      software.assigned_scores.push_back(refs.score_ref);
      software.assigned_scores.push_back(refs.pep_ref);
      rescoring_scores = refs;
    }
    IdentificationData::ProcessingSoftwareRef software_ref =
      id_data.registerProcessingSoftware(software);
    // @TODO: add suitable data processing action
    IdentificationData::ProcessingStep step(software_ref, {file_ref});
    if (route_generation)
    {
      // OMS stores ProcessingStep metadata, but does not persist DBSearchParam
      // metadata. Keep the cache contract on the processing step as well.
      step.setMetaValue("digestion:placement_model", canonical_cleavage ? "route_canonical" : "route_deferred");
      step.setMetaValue("digestion:route_max_modifications", static_cast<unsigned int>(max_variable_mods_per_oligo));
    }
    step.actions.insert(DataProcessing::IDENTIFICATION);
    IdentificationData::ProcessingStepRef step_ref;

    // get digested sequences:
    std::string decoy_pattern = getStringOption_("fdr:decoy_pattern");
    if (rescoring_enabled && decoy_pattern.empty())
    {
      OPENMS_LOG_ERROR << "Error: 'rescoring:enable' requires decoy sequences ('fdr:decoy_pattern')." << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (in_digest.empty()) // new digestion
    {
      search_param.database = in_db;
      search_param.missed_cleavages = getIntOption_("oligo:missed_cleavages");
      search_param.digestion_enzyme = configured_enzyme;
      search_param.min_length = getIntOption_("oligo:min_size");
      search_param.max_length = getIntOption_("oligo:max_size");
      if (route_generation)
      {
        search_param.setMetaValue("digestion:placement_model", canonical_cleavage ? "route_canonical" : "route_deferred");
        search_param.setMetaValue("digestion:route_max_modifications", static_cast<unsigned int>(max_variable_mods_per_oligo));
      }
      IdentificationData::SearchParamRef search_ref =
        id_data.registerDBSearchParam(search_param);
      step_ref = id_data.registerProcessingStep(step, search_ref);
      // reference this step in all following ID data items, if applicable:
      id_data.setCurrentProcessingStep(step_ref);

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
      IdentificationDataConverter::importSequences(
        id_data, fasta_db, IdentificationData::MoleculeType::RNA, decoy_pattern);
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
      for (auto parent_ref = id_data.getParentSequences().begin();
           parent_ref != id_data.getParentSequences().end(); ++parent_ref)
      {
        if (parent_ref->molecule_type != IdentificationData::MoleculeType::RNA)
        {
          continue;
        }
        NASequence parent = NASequence::fromString(parent_ref->sequence);
        ModifiedNASequenceGenerator::applyFixedModifications(internal_fixed, parent);
        vector<RNaseDigestion::DigestionProduct> products;
        if (route_generation)
        {
          // Store nucleotide sequences once. Modification-dependent products
          // are validated against their full parent after focused placement.
          const auto* rnase = dynamic_cast<const DigestionEnzymeRNA*>(search_param.digestion_enzyme);
          if (! rnase)
          {
            OPENMS_LOG_ERROR << "graph_routes requires an RNase for digestion.\n";
            return ILLEGAL_PARAMETERS;
          }
          StringList after, before;
          StringUtils::split(rnase->getCutsAfterRegEx(), ',', after);
          StringUtils::split(rnase->getCutsBeforeRegEx(), ',', before);
          const Size width = std::max<Size>(1, after.size() + before.size());
          const Size extra_cuts = std::min(parent.size(), max_variable_mods_per_oligo) * width;
          digestor.setMissedCleavages(search_param.missed_cleavages + (canonical_cleavage ? 0 : extra_cuts));
          digestor.digest(parent, products, min_oligo_length, max_oligo_length);
          for (auto& product : products)
          {
            if (product.position.first != 0 || ! parent.hasFivePrimeMod()) { product.fragment.setFivePrimeMod(nullptr); }
            if (product.position.first + product.position.second != parent.size() || ! parent.hasThreePrimeMod())
            {
              product.fragment.setThreePrimeMod(nullptr);
            }
            ModifiedNASequenceGenerator::applyFixedModifications(fixed_modifications, product.fragment);
          }
        }
        else
        {
          digestor.digestWithCleavageSensitiveMods(parent, cleavage_sensitive_groups, max_variable_mods_per_oligo, products, min_oligo_length,
                                                   max_oligo_length);
        }
        for (const auto& product : products)
        {
          registerDigestedOligo_(id_data, parent_ref, parent, product);
        }
      }

      std::string digest_out = getStringOption_("digest_out");
      if (!digest_out.empty())
      {
        OMSFile(log_type_).store(digest_out, id_data);
      }
    }
    else // load digestion results from a previous run
    {
      OPENMS_LOG_INFO << "Loading pre-digested sequence data..." << endl;
      OMSFile(log_type_).load(in_digest, id_data);
      if (route_generation && ! id_data.getObservationMatches().empty())
      {
        OPENMS_LOG_ERROR << "graph_routes requires a sequence digest from 'digest_out', not search results containing observation matches.\n";
        return ILLEGAL_PARAMETERS;
      }
      if (id_data.getDBSearchParams().empty())
      {
        OPENMS_LOG_ERROR << "No search parameter information found in 'digest' input" << endl;
        return ILLEGAL_PARAMETERS;
      }

      IdentificationData::SearchParamRef search_ref =
        id_data.getDBSearchParams().begin();
      const std::string required_digest = route_generation ? (canonical_cleavage ? "route_canonical" : "route_deferred") : "preexpanded";
      bool matching_route_digest = false;
      bool has_route_digest = false;
      Size cached_route_budget = 0;
      for (const auto& cached_step : id_data.getProcessingSteps())
      {
        if (cached_step.metaValueExists("digestion:placement_model"))
        {
          const auto model = cached_step.getMetaValue("digestion:placement_model").toString();
          has_route_digest |= model == "route_canonical" || model == "route_deferred";
        }
        if (cached_step.metaValueExists("digestion:placement_model") && cached_step.getMetaValue("digestion:placement_model") == required_digest)
        {
          matching_route_digest = true;
          if (cached_step.metaValueExists("digestion:route_max_modifications"))
          {
            cached_route_budget
              = std::max(cached_route_budget, Size(static_cast<unsigned int>(cached_step.getMetaValue("digestion:route_max_modifications"))));
          }
        }
      }
      if (! route_generation && has_route_digest)
      {
        OPENMS_LOG_ERROR << "Enumeration/tag modes require a pre-expanded digest; create a fresh digest from the FASTA.\n";
        return ILLEGAL_PARAMETERS;
      }
      if (route_generation && ! matching_route_digest)
      {
        OPENMS_LOG_ERROR << "graph_routes requires a matching deferred/canonical digest; create a fresh digest from the FASTA.\n";
        return ILLEGAL_PARAMETERS;
      }
      if (route_generation
          && (search_ref->fixed_mods != search_param.fixed_mods || (! canonical_cleavage && cached_route_budget < max_variable_mods_per_oligo)))
      {
        OPENMS_LOG_ERROR << "Deferred digest fixed modifications/budget do not cover this search; create a fresh digest.\n";
        return ILLEGAL_PARAMETERS;
      }
      search_param.digestion_enzyme = search_ref->digestion_enzyme;
      // Route validation re-digests products with the cached budget.
      search_param.missed_cleavages = search_ref->missed_cleavages;
      cleavage_sensitive_groups = {};
      cleavage_sensitive_modifications.clear();
      cleavage_sensitive_mod_codes.clear();
      post_digest_variable_modifications = variable_modifications;
      if (search_ref->digestion_enzyme)
      {
        RNaseDigestion cleavage_sensitive_digestor;
        cleavage_sensitive_digestor.setEnzyme(search_ref->digestion_enzyme);
        cleavage_sensitive_groups =
          cleavage_sensitive_digestor.inferCleavageSensitiveMods(
            variable_modifications);
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
      }

      // Preserve digestion/database provenance, but record the parameters used
      // by this search rather than the tolerances of the cached digest producer.
      IdentificationData::DBSearchParam current_search = search_param;
      current_search.database = search_ref->database;
      current_search.database_version = search_ref->database_version;
      current_search.taxonomy = search_ref->taxonomy;
      current_search.digestion_enzyme = search_ref->digestion_enzyme;
      current_search.enzyme_term_specificity = search_ref->enzyme_term_specificity;
      current_search.missed_cleavages = search_ref->missed_cleavages;
      current_search.min_length = search_ref->min_length;
      current_search.max_length = search_ref->max_length;
      std::vector<std::string> cached_keys;
      search_ref->getKeys(cached_keys);
      for (const auto& key : cached_keys)
      {
        if (key.starts_with("digestion:")) current_search.setMetaValue(key, search_ref->getMetaValue(key));
      }
      search_ref = id_data.registerDBSearchParam(current_search);
      // Fresh and cached searches both record the identification action.
      step_ref = id_data.registerProcessingStep(step, search_ref);
      // reference this step in all following ID data items:
      id_data.setCurrentProcessingStep(step_ref);
    }
    Size n_nucleic_acids = id_data.getParentSequences().size();

    if (!decoy_pattern.empty())
    {
      bool no_decoys = none_of(
        id_data.getParentSequences().begin(),
        id_data.getParentSequences().end(),
        [](const IdentificationData::ParentSequence& p){ return p.is_decoy; });
      if (no_decoys)
      {
        OPENMS_LOG_ERROR << "Error: 'fdr:decoy_pattern' is set, but no decoy sequences were found" << endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    progresslogger.startProgress(0, 1, "filtering spectra...");
    // @TODO: move this into the loop below (run only when checks pass)
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
    const bool native_labile_scaffolds = getFlag_("scoring:labile_scaffolds");
    if (native_labile_scaffolds && (! route_generation || scoring_method == "route"))
    {
      OPENMS_LOG_ERROR << "scoring:labile_scaffolds requires graph_routes with MVH or HyperScore.\n";
      return ILLEGAL_PARAMETERS;
    }
    const bool diagnostic_filter = getFlag_("diagnostic:filter");
    const bool diagnostic_evaluate = diagnostic_filter || route_generation || scoring_method == "route";
    const bool explicit_diagnostic_filter = getStringOption_("diagnostic:filter_scope") == "explicit";
    std::set<std::string> explicit_diagnostic_codes;
    std::vector<RNAModificationDiagnostics::Evidence> diagnostic_evidence;
    if (diagnostic_evaluate)
    {
      std::map<std::string, std::vector<double>> fixed_diagnostics;
      for (const auto& item : getStringList_("diagnostic:fixed_ions"))
      {
        const auto separator = item.find('=');
        if (separator == std::string::npos || separator == 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Diagnostic rules require code=mz");
        }
        explicit_diagnostic_codes.insert(item.substr(0, separator));
        fixed_diagnostics[item.substr(0, separator)].push_back(StringUtils::toDouble(item.substr(separator + 1)));
      }
      RNAModificationDiagnostics diagnostics(variable_modifications, fixed_diagnostics, negative_mode, search_param.fragment_mass_tolerance,
                                             search_param.fragment_tolerance_ppm, getDoubleOption_("diagnostic:min_relative_intensity"),
                                             getIntOption_("diagnostic:max_charge"));
      // Sort once before attaching evidence; preprocessing keeps spectrum order.
      spectra.sortSpectra(false);
      for (auto& spectrum : spectra)
      {
        auto evidence = diagnostics.evaluate(spectrum);
        spectrum.setMetaValue("diagnostic_allowed_mods", StringList(evidence.allowed.begin(), evidence.allowed.end()));
        spectrum.setMetaValue("diagnostic_labile_mods", StringList(evidence.labile.begin(), evidence.labile.end()));
        StringList diagnostic_matches;
        for (const auto& [code, mz] : evidence.matched_mz)
        {
          diagnostic_matches.push_back(code + "=" + StringUtils::toStr(mz));
        }
        spectrum.setMetaValue("diagnostic_matches", diagnostic_matches);
        StringList diagnostic_charges, diagnostic_isotopes;
        for (const auto& [code, charge] : evidence.matched_charge)
        {
          diagnostic_charges.push_back(code + "=" + StringUtils::toStr(charge));
        }
        for (const auto& [code, mz] : evidence.isotope_mz)
        {
          diagnostic_isotopes.push_back(code + "=" + StringUtils::toStr(mz));
        }
        spectrum.setMetaValue("diagnostic_match_charges", diagnostic_charges);
        spectrum.setMetaValue("diagnostic_isotope_matches", diagnostic_isotopes);
        diagnostic_evidence.push_back(std::move(evidence));
      }
    }
    preprocessSpectra_(spectra, search_param.fragment_mass_tolerance, search_param.fragment_tolerance_ppm, single_charge_spectra, negative_mode,
                       min_charge, max_charge, include_unknown_charge, use_window_mower, window_size, window_peakcount, window_movetype, use_nlargest,
                       nlargest_n, remove_precursor, precursor_mass_tolerance, precursor_tolerance_ppm, precursor_peak_isotopes, deisotoping,
                       flash_min_cos, flash_min_snr, route_generation || (scoring_method == "route"),
                       getIntOption_("preprocessing:legacy:min_isopeaks"), getDoubleOption_("preprocessing:legacy:mass_tolerance"),
                       getFlag_("preprocessing:legacy:sum_intensities"), getIntOption_("preprocessing:legacy:unknown_max_charge"));
    if (diagnostic_evaluate)
    {
      // Recover evidence from spectrum-owned metadata after preprocessing sorts;
      // equal-RT spectra and duplicate native IDs must not exchange permissions.
      diagnostic_evidence.clear();
      for (const auto& spectrum : spectra)
      {
        RNAModificationDiagnostics::Evidence evidence;
        for (const auto* residue : variable_modifications)
        {
          if (diagnostic_filter && (! explicit_diagnostic_filter || explicit_diagnostic_codes.contains(residue->getCode())))
          {
            evidence.tested.insert(residue->getCode());
          }
        }
        for (const auto& code : spectrum.getMetaValue("diagnostic_allowed_mods").toStringList())
        {
          evidence.allowed.insert(code);
        }
        for (const auto& code : spectrum.getMetaValue("diagnostic_labile_mods").toStringList())
        {
          evidence.labile.insert(code);
        }
        diagnostic_evidence.push_back(std::move(evidence));
      }
    }
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
      const bool flash_matched = correct_precursor && s_it->metaValueExists("flash_precursor_status") &&
        (s_it->getMetaValue("flash_precursor_status") == "corrected" ||
         s_it->getMetaValue("flash_precursor_status") == "monoisotopic");
      const IntList scan_isotopes = flash_matched ? IntList{0} : precursor_isotopes;


      for (Int precursor_charge : possible_charges)
      {
        precursor_charge = abs(precursor_charge); // adjust for neg. mode

        // calculate precursor mass (optionally corrected for adducts and peak
        // misassignment) and map it to MS scan index etc.:
        for (const auto& adduct_pair : adduct_masses)
        {
          for (Int isotope_number : scan_isotopes)
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

    vector<RNASpectrumGraph::Result> spectrum_graphs;
    std::shared_ptr<const RNASpectrumGraph::Alphabet> graph_alphabet;
    if (graph_enabled)
    {
      RNASpectrumGraph::Options graph_options;
      graph_options.tolerance = search_param.fragment_mass_tolerance;
      graph_options.tolerance_ppm = search_param.fragment_tolerance_ppm;
      graph_options.min_tag_length = getIntOption_("graph:min_tag_length");
      graph_options.max_tag_length = getIntOption_("graph:max_tag_length");
      graph_options.max_step = getIntOption_("graph:max_step");
      graph_options.max_peaks = getIntOption_("graph:max_peaks");
      graph_options.max_charge = getIntOption_("graph:max_charge");
      graph_options.max_edges = getIntOption_("graph:max_edges");
      graph_options.max_query_edges = getIntOption_("graph:max_query_edges");
      graph_options.max_query_states = getIntOption_("graph:max_query_states");
      const Size total_edge_budget = getIntOption_("graph:max_total_edges");
      set<ConstRibonucleotidePtr> alphabet = variable_modifications;
      alphabet.insert(fixed_modifications.begin(), fixed_modifications.end());
      for (const auto& code : {"A", "C", "G", "U"})
      {
        alphabet.insert(RibonucleotideDB::getInstance()->getRibonucleotide(code));
      }
      const RNASpectrumGraph graph_builder(alphabet, graph_options);
      spectrum_graphs.resize(spectra.size());
      Size total_edges = 0;
      map<std::string, Size> fallbacks;
      progresslogger.startProgress(0, spectra.size(), "building RNA mass-tag graphs...");
      // Deterministic spectrum order makes the shared cache budget independent
      // of thread scheduling and applies the same graph to targets and decoys.
      for (Size scan = 0; scan < spectra.size(); ++scan)
      {
        progresslogger.setProgress(scan);
        auto& graph = spectrum_graphs[scan];
        if (spectra[scan].getMSLevel() != 2 || spectra[scan].getPrecursors().size() != 1) { graph.fallback_reason = "invalid_ms2_precursor"; }
        else
        {
          Size charge = std::abs(spectra[scan].getPrecursors()[0].getCharge());
          if (charge == 0)
          {
            charge = std::max(std::abs(min_charge), std::abs(max_charge));
          }
          graph = graph_builder.build(spectra[scan], single_charge_spectra ? 1 : charge);
          if (graph.fallback_reason.empty())
          {
            // Enumerating only the remaining spectra would make candidate
            // selection, and therefore FDR, depend on the spectrum order.
            if (graph.edges > total_edge_budget - total_edges)
            {
              progresslogger.endProgress();
              OPENMS_LOG_ERROR << "RNA graphs exceed 'graph:max_total_edges' (" << total_edge_budget << ") at spectrum " << (scan + 1) << " of "
                               << spectra.size() << ". Raise the bound, lower 'graph:max_peaks' or use 'candidates:method enumeration'.\n";
              return ILLEGAL_PARAMETERS;
            }
            total_edges += graph.edges;
          }
        }
        if (! graph.fallback_reason.empty())
        {
          ++fallbacks[graph.fallback_reason];
        }
        if (! graph_alphabet && graph.alphabet) { graph_alphabet = graph.alphabet; }
      }
      progresslogger.endProgress();
      OPENMS_LOG_INFO << "RNA graph cached edges: " << total_edges << '\n';
      for (const auto& [reason, count] : fallbacks)
      {
        OPENMS_LOG_INFO << "RNA graph fallback " << reason << ": " << count << '\n';
      }
      if (! graph_out.empty())
      {
        SVOutStream output(graph_out, "\t", "_");
        output << "spectrum" << "peaks" << "vertices" << "edges" << "alphabet" << "fallback" << endl;
        for (Size scan = 0; scan < spectra.size(); ++scan)
        {
          const auto& graph = spectrum_graphs[scan];
          output << spectra[scan].getNativeID() << graph.peaks << graph.vertices << graph.edges
                 << (graph.alphabet ? graph.alphabet->code_ids.size() : 0) << graph.fallback_reason << endl;
        }
      }
    }

    vector<HitsByScore> annotated_hits(spectra.size());
    MSExperiment exp_ms2_spectra, theo_ms2_spectra; // debug output

    std::string msg = "scoring oligonucleotide models against spectra...";
    progresslogger.startProgress(0, id_data.getIdentifiedOligos().size(), msg);
    Size hit_counter = 0;
    Size graph_prefiltered = 0, graph_rejected = 0, graph_supported = 0, graph_fallback = 0;
    Size graph_verified = 0, graph_mismatches = 0;
    Size placement_nodes = 0, placement_pruned = 0, placement_emitted = 0, placement_fallbacks = 0;
    Size placement_reference = 0, placement_expected = 0, placement_missing = 0, placement_extra = 0;

    // Deferred products are search scaffolds, not identified output molecules.
    // Preserve their references and the original parent references in a separate
    // container. Otherwise registration can merge a filtered candidate back into
    // its deferred scaffold and reintroduce terminally incompatible parents.
    IdentificationData::IdentifiedOligos deferred_digest;
    if (route_generation)
    {
      deferred_digest = id_data.getIdentifiedOligos();
      id_data.cleanup(true, false); // remove unmatched scaffolds, keep parent references
    }
    const auto& search_oligos = route_generation ? deferred_digest : id_data.getIdentifiedOligos();
    vector<IdentificationData::IdentifiedOligoRef> digest;
    digest.reserve(search_oligos.size());
    for (IdentificationData::IdentifiedOligoRef it = search_oligos.begin(); it != search_oligos.end(); ++it)
    {
      digest.push_back(it);
    }

    Int base_charge = negative_mode ? -1 : 1;
    const auto blacklist_mz = getDoubleList_("scoring:blacklist_mz");
    const double blacklist_tolerance = getDoubleOption_("scoring:blacklist_tolerance");
    const Size num_intensity_classes = getIntOption_("scoring:num_intensity_classes");
    const double tic_fraction = getDoubleOption_("scoring:tic_fraction");
    const bool use_mass_accuracy = getFlag_("scoring:use_mass_accuracy");
    const bool route_base_loss = getFlag_("scoring:route_base_loss");
    const bool single_peak_credit = getStringOption_("scoring:single_peak_credit") == "true";
    const double min_coverage = getDoubleOption_("report:require_coverage");
    const Size route_max_charge = getIntOption_("graph:max_charge");
    RNARouteScorer::Options route_score_options;
    route_score_options.tolerance = search_param.fragment_mass_tolerance;
    route_score_options.tolerance_ppm = search_param.fragment_tolerance_ppm;
    route_score_options.local_radius = getDoubleOption_("scoring:route_local_radius");
    route_score_options.mass_accuracy_tolerance = getDoubleOption_("scoring:route_mass_accuracy_tolerance");
    route_score_options.local_peaks = getIntOption_("scoring:route_local_peaks");
    route_score_options.min_nodes = getIntOption_("scoring:route_min_nodes");
    route_score_options.max_unknown_charge = single_charge_spectra ? 1 : route_max_charge;
    route_score_options.max_states = getIntOption_("routes:max_states");
    vector<MSSpectrum> route_scoring_spectra;
    vector<RNARouteScorer::PreparedSpectrum> route_scoring_contexts;
    if (scoring_method == "route")
    {
      route_scoring_spectra.assign(spectra.begin(), spectra.end());
      const Size max_peaks = getIntOption_(route_generation ? "routes:max_peaks" : "graph:max_peaks");
      for (auto& spectrum : route_scoring_spectra)
      {
        spectrum.select(RNASpectrumGraph::selectPeakIndices(spectrum, max_peaks));
        route_scoring_contexts.push_back(RNARouteScorer::prepareSpectrum(spectrum, negative_mode, route_score_options));
      }
    }

    std::atomic<bool> route_score_limit_exceeded {false};
    std::string route_score_failure; // first failing pair; guarded by an OpenMP critical section
    if (route_generation)
    {
      vector<NASequence> route_database;
      for (const auto& reference : digest)
      {
        route_database.push_back(reference->sequence);
      }
      RNARouteCandidateGenerator::Options options;
      options.tolerance = search_param.fragment_mass_tolerance;
      options.tolerance_ppm = search_param.fragment_tolerance_ppm;
      options.precursor_tolerance = search_param.precursor_mass_tolerance;
      options.precursor_tolerance_ppm = search_param.precursor_tolerance_ppm;
      options.max_modifications = max_variable_mods_per_oligo;
      options.max_length = getIntOption_("routes:max_length");
      options.max_peaks = getIntOption_("routes:max_peaks");
      options.max_states = getIntOption_("routes:max_states");
      options.max_completion_nodes = getIntOption_("routes:max_completion_nodes");
      options.max_candidates = getIntOption_("routes:max_candidates");
      options.max_supports = getIntOption_("routes:max_supports");
      options.min_route_nodes = getIntOption_("routes:min_nodes");
      const auto completion = getStringOption_("routes:completion");
      if ((completion == "labile" || completion == "bounded_labile") && ! variable_modifications.empty()
          && getStringList_("diagnostic:fixed_ions").empty())
      {
        OPENMS_LOG_WARN << "Warning: 'routes:completion " << completion << "' completes uncovered positions only with labile modifications, "
                        << "and no 'diagnostic:fixed_ions' rule declares one. Variable modifications without an intact graph placement "
                        << "are not searched; set 'diagnostic:fixed_ions' or use 'routes:completion bounded'.\n";
      }
      options.complete_unanchored_modifications = completion == "bounded" || completion == "all";
      if (completion == "bounded" || completion == "bounded_labile")
      {
        options.max_unanchored_modifications = getIntOption_("routes:max_unanchored_modifications");
      }
      options.index_precursor_mass = getStringOption_("routes:precursor_index") == "true";
      options.allow_full_length_routes = getStringOption_("routes:allow_full_length") == "true";
      // Terminal groups follow the RNase and the parent position of each
      // product, as in enumeration. One sequence can occur at several
      // positions: the generator receives their union, and each parent match
      // is checked again once a candidate carries concrete termini.
      const bool enzyme_termini = getStringOption_("routes:terminal_states") == "enzyme";
      const auto* route_rnase = dynamic_cast<const DigestionEnzymeRNA*>(search_param.digestion_enzyme);
      if ((enzyme_termini || ! canonical_cleavage) && ! route_rnase)
      {
        OPENMS_LOG_ERROR << "graph_routes requires an RNase, unless 'routes:terminal_states' is 'all' and 'oligo:cleavage_model' is 'canonical'.\n";
        return ILLEGAL_PARAMETERS;
      }
      Size five_gain = 0, three_gain = 0;
      if (enzyme_termini)
      {
        const std::string five_code = route_rnase->getFivePrimeGain(), three_code = route_rnase->getThreePrimeGain();
        if ((! five_code.empty() && five_code != "p") || (! three_code.empty() && three_code != "p" && three_code != "c"))
        {
          OPENMS_LOG_ERROR << "graph_routes supports OH, phosphate and 3' cyclic phosphate termini; the terminal gains of '"
                           << route_rnase->getName() << "' require enumeration.\n";
          return ILLEGAL_PARAMETERS;
        }
        five_gain = five_code.empty() ? 0 : 1;
        three_gain = three_code.empty() ? 0 : (three_code == "p" ? 1 : 2);
      }
      // Index 3*five+three, in the order used by RNARouteCandidateGenerator.
      const auto ownTerminalState = [](const NASequence& sequence) {
        const Size three = ! sequence.hasThreePrimeMod() ? 0 : (sequence.getThreePrimeMod()->getCode() == "3'-p" ? 1 : 2);
        return 3 * Size(sequence.hasFivePrimeMod()) + three;
      };
      // As in enumeration, an end takes the group gained from the RNase. An
      // end without one (a parent end, or no gain) keeps the group stored
      // with the digest: the parent's own or a fixed terminal modification.
      const auto enzymeTerminalState = [&](const NASequence& stored, const IdentificationData::ParentMatch& match, Size parent_size) {
        const Size own = ownTerminalState(stored);
        return 3 * (match.start_pos == 0 || five_gain == 0 ? own / 3 : five_gain)
               + (match.end_pos + 1 >= parent_size || three_gain == 0 ? own % 3 : three_gain);
      };
      // Parse each parent once; products are validated against it below.
      std::map<IdentificationData::ParentSequenceRef, NASequence> route_parents;
      vector<Byte> terminal_masks;
      if (enzyme_termini || ! canonical_cleavage)
      {
        if (enzyme_termini) { terminal_masks.assign(digest.size(), 0); }
        for (Size index = 0; index < digest.size(); ++index)
        {
          for (const auto& [parent_ref, matches] : digest[index]->parent_matches)
          {
            auto found = route_parents.find(parent_ref);
            if (found == route_parents.end())
            {
              NASequence parent = NASequence::fromString(parent_ref->sequence);
              ModifiedNASequenceGenerator::applyFixedModifications(fixed_modifications, parent);
              found = route_parents.emplace(parent_ref, std::move(parent)).first;
            }
            if (! enzyme_termini) { continue; }
            for (const auto& match : matches)
            {
              terminal_masks[index] |= Byte(1) << enzymeTerminalState(digest[index]->sequence, match, found->second.size());
            }
          }
        }
      }
      // A cleavage site depends only on the residues matched by the enzyme
      // patterns around it. Digest the product with that much parent context
      // (at least one residue, so a window edge is never mistaken for a cut).
      Size context_before = 1, context_after = 1;
      if (route_rnase)
      {
        StringList after, before;
        StringUtils::split(route_rnase->getCutsAfterRegEx(), ',', after);
        StringUtils::split(route_rnase->getCutsBeforeRegEx(), ',', before);
        context_before = std::max<Size>(1, after.size());
        context_after = std::max<Size>(1, before.size());
      }
      const auto cleavageCompatible = [&](const RNaseDigestion& validator, const NASequence& parent, Size start, const NASequence& candidate) {
        const Size end = start + candidate.size();
        if (end > parent.size()) { return false; }
        const Size first = start - std::min(start, context_before), last = end + std::min(parent.size() - end, context_after);
        std::vector<const Ribonucleotide*> residues;
        residues.reserve(last - first);
        for (Size pos = first; pos < last; ++pos)
        {
          residues.push_back(pos >= start && pos < end ? candidate[pos - start] : parent[pos]);
        }
        const NASequence window(std::move(residues), nullptr, nullptr);
        const auto positions = validator.getFragmentPositions(window, candidate.size(), candidate.size());
        return std::find(positions.begin(), positions.end(), std::pair<Size, Size> {start - first, candidate.size()}) != positions.end();
      };
      const RNARouteCandidateGenerator generator(route_database, variable_modifications, options, terminal_masks);
      std::unique_ptr<RNARouteCandidateGenerator> route_reference;
      if (getFlag_("routes:verify_candidates"))
      {
        auto reference_options = options;
        reference_options.index_precursor_mass = false;
        reference_options.compress_route_states = false;
        route_reference = std::make_unique<RNARouteCandidateGenerator>(route_database, variable_modifications, reference_options, terminal_masks);
      }
      using PrecursorEntry = std::pair<const double, PrecursorInfo>;
      vector<vector<const PrecursorEntry*>> by_scan(spectra.size());
      for (const auto& precursor : precursor_mass_map)
      {
        by_scan[precursor.second.scan_index].push_back(&precursor);
      }
      vector<RNARouteCandidateGenerator::Result> route_results(spectra.size());
      vector<std::string> failures(spectra.size());
      std::atomic<bool> route_incomplete {false};
#pragma omp parallel for schedule(dynamic) reduction(+ : hit_counter)
      for (SignedSize scan = 0; scan < SignedSize(spectra.size()); ++scan)
      {
        if (by_scan[scan].empty()) { continue; }
        if (route_incomplete.load(std::memory_order_relaxed))
        {
          failures[scan] = "not_attempted_after_incomplete_search";
          route_results[scan].complete = false;
          route_results[scan].reason = failures[scan];
          continue;
        }
        try
        {
          vector<double> masses;
          Size max_fragment_charge = 1;
          for (const auto* precursor : by_scan[scan])
          {
            masses.push_back(precursor->first);
            max_fragment_charge = std::max(max_fragment_charge, Size(precursor->second.charge));
          }
          const auto& evidence = diagnostic_evidence[scan];
          auto generated = generator.generate(spectra[scan], masses, evidence, negative_mode,
                                              single_charge_spectra ? 1 : std::min(max_fragment_charge, route_max_charge), &route_incomplete);
          if (generated.complete && route_reference)
          {
            const auto reference
              = route_reference->generate(spectra[scan], masses, evidence, negative_mode,
                                          single_charge_spectra ? 1 : std::min(max_fragment_charge, route_max_charge), &route_incomplete);
            const auto keys = [](const auto& candidates) {
              std::vector<std::tuple<Size, Size, std::string, std::vector<std::string>>> result;
              for (const auto& candidate : candidates)
              {
                std::vector<std::string> scaffolds;
                for (const auto& scaffold : candidate.support_scaffolds)
                {
                  scaffolds.push_back(scaffold->label);
                }
                std::sort(scaffolds.begin(), scaffolds.end());
                result.emplace_back(candidate.database_index, candidate.precursor_index, candidate.sequence.toString(), std::move(scaffolds));
              }
              std::sort(result.begin(), result.end());
              return result;
            };
            if (! reference.complete || keys(reference.candidates) != keys(generated.candidates))
            {
              generated.complete = false;
              generated.reason = reference.complete ? "candidate_multiset_mismatch" : "reference:" + reference.reason;
            }
          }
          if (! generated.complete)
          {
            failures[scan] = generated.reason;
            route_incomplete.store(true, std::memory_order_relaxed);
          }
          else
          {
            std::map<std::pair<const RNARouteCandidateGenerator::Scaffold*, Size>, RNARouteScorer::PreparedRoutes> support_routes;
            struct NativeModel
            {
              double score;
              std::vector<PeptideHit::PeakAnnotation> annotations;
            };
            std::map<std::pair<const RNARouteCandidateGenerator::Scaffold*, Size>, NativeModel> native_support_models;
            RNaseDigestion validator;
            if (! canonical_cleavage)
            {
              validator.setEnzyme(search_param.digestion_enzyme);
              validator.setMissedCleavages(search_param.missed_cleavages);
            }
            const auto score_native = [&](const MSSpectrum& model) {
              NativeModel result;
              result.score
                = scoreNativeSpectrum_(scoring_method, search_param.fragment_mass_tolerance, search_param.fragment_tolerance_ppm,
                                                                           spectra[scan], model, result.annotations, num_intensity_classes,
                                                                           tic_fraction, use_mass_accuracy, single_peak_credit);
              return result;
            };
            for (const auto& placed : generated.candidates)
            {
              if (route_incomplete.load(std::memory_order_relaxed))
              {
                failures[scan] = "cancelled_after_incomplete_search";
                break;
              }
              const auto* precursor = by_scan[scan][placed.precursor_index];
              const auto& candidate = placed.sequence;
              std::optional<IdentificationData::IdentifiedOligo> parent_override;
              if (enzyme_termini || ! canonical_cleavage)
              {
                parent_override = *digest[placed.database_index];
                parent_override->parent_matches.clear();
                const Size terminal_state = ownTerminalState(candidate);
                for (const auto& [parent_ref, matches] : digest[placed.database_index]->parent_matches)
                {
                  const NASequence& parent = route_parents.at(parent_ref);
                  for (const auto& match : matches)
                  {
                    if (enzyme_termini && enzymeTerminalState(digest[placed.database_index]->sequence, match, parent.size()) != terminal_state)
                    {
                      continue;
                    }
                    if (canonical_cleavage || cleavageCompatible(validator, parent, match.start_pos, candidate))
                    {
                      parent_override->parent_matches[parent_ref].insert(match);
                    }
                  }
                }
                if (parent_override->parent_matches.empty()) { continue; }
              }
              MSSpectrum model;
              spectrum_generator.getSpectrum(model, candidate, base_charge, precursor->second.charge * base_charge);
              filterBlacklistedIons_(model, blacklist_mz, blacklist_tolerance);
              vector<PeptideHit::PeakAnnotation> annotations;
              RNARouteScorer::Result route_evidence;
              std::string native_scoring_scaffold;
              double score;
              if (scoring_method == "route")
              {
                const auto query = RNATagCandidateFilter::prepare(candidate, nullptr, model, selected_ions, route_base_loss);
                const auto intact = RNARouteScorer::match(query, route_scoring_spectra[scan], route_scoring_contexts[scan], negative_mode,
                                                          route_score_options, candidate.toString());
                std::vector<const RNARouteScorer::PreparedRoutes*> routes {&intact};
                for (const auto& scaffold : placed.support_scaffolds)
                {
                  const auto key = std::pair(scaffold.get(), precursor->second.charge);
                  auto found = support_routes.find(key);
                  if (found == support_routes.end())
                  {
                    MSSpectrum support_model;
                    spectrum_generator.getSpectrum(support_model, scaffold->sequence, base_charge, precursor->second.charge * base_charge);
                    filterBlacklistedIons_(support_model, blacklist_mz, blacklist_tolerance);
                    const auto support_query
                      = RNATagCandidateFilter::prepare(scaffold->sequence, nullptr, support_model, selected_ions, route_base_loss);
                    auto matched = RNARouteScorer::match(support_query, route_scoring_spectra[scan], route_scoring_contexts[scan], negative_mode,
                                                         route_score_options, scaffold->label);
                    found = support_routes.emplace(key, std::move(matched)).first;
                  }
                  routes.push_back(&found->second);
                }
                route_evidence = RNARouteScorer::score(routes, route_scoring_spectra[scan], route_scoring_contexts[scan], negative_mode);
                score = route_evidence.score;
                annotations = route_evidence.annotations;
              }
              else
              {
                auto best = score_native(model);
                if (native_labile_scaffolds)
                {
                  native_scoring_scaffold = candidate.toString();
                  for (const auto& scaffold : placed.support_scaffolds)
                  {
                    const auto key = std::pair(scaffold.get(), precursor->second.charge);
                    auto found = native_support_models.find(key);
                    if (found == native_support_models.end())
                    {
                      MSSpectrum support_model;
                      spectrum_generator.getSpectrum(support_model, scaffold->sequence, base_charge, precursor->second.charge * base_charge);
                      filterBlacklistedIons_(support_model, blacklist_mz, blacklist_tolerance);
                      found = native_support_models.emplace(key, score_native(support_model)).first;
                    }
                    // Keep one coherent model. Do not pool ions from mutually
                    // incompatible scaffold interpretations.
                    if (std::isfinite(found->second.score) && (! std::isfinite(best.score) || found->second.score > best.score))
                    {
                      best = found->second;
                      native_scoring_scaffold = scaffold->label;
                    }
                  }
                }
                score = best.score;
                annotations = std::move(best.annotations);
              }
              if (! std::isfinite(score) || score < 1e-16 || (min_coverage > 0 && calculateFragmentCoverage_(candidate, annotations) < min_coverage))
              {
                continue;
              }
              ++hit_counter;
              AnnotatedHit hit;
              hit.oligo_ref = digest[placed.database_index];
              hit.sequence = candidate;
              hit.annotations = std::move(annotations);
              hit.precursor_error_ppm = (precursor->first - candidate.getMonoWeight()) / candidate.getMonoWeight() * 1e6;
              hit.precursor_ref = &precursor->second;
              hit.route_evidence = std::move(route_evidence);
              hit.native_scoring_scaffold = std::move(native_scoring_scaffold);
              hit.parent_override = std::move(parent_override);
              retainTopHit_(annotated_hits[scan], score, std::move(hit), search_top_hits);
            }
          }
          generated.candidates.clear();
          generated.candidates.shrink_to_fit();
          route_results[scan] = std::move(generated);
        }
        catch (const std::exception& error)
        {
          failures[scan] = error.what();
          route_incomplete.store(true, std::memory_order_relaxed);
        }
      }
      std::ofstream diagnostics;
      if (! graph_out.empty())
      {
        diagnostics.open(graph_out);
        diagnostics << "scan_index\tnative_id\tgraph_nodes\tgraph_edges\troute_states\troutes\tcompletion_nodes\tsupport_variants\tstatus\n";
      }
      bool complete = true;
      Size cancelled_scans = 0;
      for (Size scan = 0; scan < spectra.size(); ++scan)
      {
        if (by_scan[scan].empty()) { continue; }
        const auto& result = route_results[scan];
        if (diagnostics)
        {
          diagnostics << scan << '\t' << spectra[scan].getNativeID() << '\t' << result.graph_nodes << '\t' << result.graph_edges << '\t'
                      << result.route_states << '\t' << result.routes << '\t' << result.completion_nodes << '\t' << result.support_variants << '\t'
                      << (failures[scan].empty() ? "complete" : failures[scan]) << '\n';
        }
        if (! failures[scan].empty())
        {
          if (failures[scan] == "not_attempted_after_incomplete_search" || failures[scan].ends_with("cancelled_after_incomplete_search"))
          {
            ++cancelled_scans;
          }
          else
          {
            // Emit root causes from the main thread: worker-local log buffers
            // can otherwise hide the reason while only cancellation is visible.
            static const std::map<std::string, std::string> limit_options {{"route_state_limit", "routes:max_states"},
                                                                            {"support_limit", "routes:max_supports"},
                                                                            {"candidate_limit", "routes:max_candidates"},
                                                                            {"completion_node_limit", "routes:max_completion_nodes"}};
            const auto limit = limit_options.find(failures[scan]);
            OPENMS_LOG_ERROR << "Incomplete route search for " << spectra[scan].getNativeID() << ": " << failures[scan]
                             << (limit == limit_options.end() ? "" : " (raise '" + limit->second + "')") << '\n';
          }
          complete = false;
        }
      }
      if (cancelled_scans) { OPENMS_LOG_INFO << "Cancelled/unattempted route spectra after incomplete search: " << cancelled_scans << '\n'; }
      if (! complete) { return UNEXPECTED_RESULT; }
    }
    else
    {

// shorter oligos take (possibly much) less time to process than longer ones;
// due to the sorting order of "NASequence", they also appear earlier in the
// container - therefore use dynamic scheduling to distribute work evenly:
#pragma omp parallel for schedule(dynamic) reduction(                                                                                          \
    + : placement_nodes, placement_pruned, placement_emitted, placement_fallbacks, placement_reference, placement_expected, placement_missing, \
      placement_extra, graph_prefiltered, graph_rejected, graph_supported, graph_fallback, graph_verified, graph_mismatches, hit_counter)
    for (SignedSize index = 0; index < SignedSize(digest.size()); ++index)
    {
      if (route_score_limit_exceeded.load(std::memory_order_relaxed)) { continue; }
      IF_MASTERTHREAD
      {
        progresslogger.setProgress(index);
      }

      IdentificationData::IdentifiedOligoRef oligo_ref = digest[index];
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
      if (graph_generate)
      {
        ModifiedNASequenceGenerator::PlacementConstraints constraints;
        constraints.use_average_mass = use_avg_mass;
        constraints.max_nodes = max_placement_nodes;
        constraints.accepts = [&](double lower_mass, double upper_mass, const auto& alternatives) {
          auto first = precursor_mass_map.begin(), past = precursor_mass_map.end();
          if (std::isfinite(lower_mass) && std::isfinite(upper_mass))
          {
            const double tolerance = search_param.precursor_mass_tolerance;
            double low = lower_mass - tolerance, high = upper_mass + tolerance;
            if (search_param.precursor_tolerance_ppm)
            {
              const double fraction = tolerance * 1e-6;
              low = std::min(lower_mass * (1.0 - fraction), upper_mass * (1.0 - fraction));
              high = std::max(lower_mass * (1.0 + fraction), upper_mass * (1.0 + fraction));
            }
            first = precursor_mass_map.lower_bound(low);
            past = precursor_mass_map.upper_bound(high);
          }
          return std::any_of(first, past, [&](const auto& precursor) {
            return RNATagCandidateFilter::isCompatible(alternatives, spectrum_graphs[precursor.second.scan_index], selected_ions,
                                                       precursor.second.charge);
          });
        };
        const auto statistics = ModifiedNASequenceGenerator::applyVariableModificationsConstrained(
          post_digest_variable_modifications, ns, remaining_variable_mods, all_modified_oligos, constraints, true);
        placement_nodes += statistics.visited_nodes;
        placement_pruned += statistics.pruned_branches;
        placement_emitted += statistics.emitted_variants;
        placement_fallbacks += statistics.enumeration_fallback;
        if (placement_audit)
        {
          vector<NASequence> reference;
          ModifiedNASequenceGenerator::applyVariableModifications(post_digest_variable_modifications, ns, remaining_variable_mods, reference, true);
          placement_reference += reference.size();
          const auto eligible = [&](const NASequence& candidate) {
            const double mass = use_avg_mass ? candidate.getAverageWeight() : candidate.getMonoWeight();
            const double tolerance = search_param.precursor_mass_tolerance * (search_param.precursor_tolerance_ppm ? mass * 1e-6 : 1.0);
            const auto first = precursor_mass_map.lower_bound(mass - tolerance), past = precursor_mass_map.upper_bound(mass + tolerance);
            return std::any_of(first, past, [&](const auto& precursor) {
              return RNATagCandidateFilter::isCompatible(candidate, spectrum_graphs[precursor.second.scan_index], selected_ions,
                                                         precursor.second.charge);
            });
          };
          vector<NASequence> expected, actual;
          std::copy_if(reference.begin(), reference.end(), std::back_inserter(expected), eligible);
          std::copy_if(all_modified_oligos.begin(), all_modified_oligos.end(), std::back_inserter(actual), eligible);
          std::sort(expected.begin(), expected.end());
          std::sort(actual.begin(), actual.end());
          placement_expected += expected.size();
          vector<NASequence> differences;
          std::set_difference(expected.begin(), expected.end(), actual.begin(), actual.end(), std::back_inserter(differences));
          placement_missing += differences.size();
          differences.clear();
          std::set_difference(actual.begin(), actual.end(), expected.begin(), expected.end(), std::back_inserter(differences));
          placement_extra += differences.size();
          // Audit scores the exhaustive reference through the phase-one filter.
          all_modified_oligos = std::move(reference);
        }
      }
      else
      {
        ModifiedNASequenceGenerator::applyVariableModifications(post_digest_variable_modifications, ns, remaining_variable_mods, all_modified_oligos,
                                                                true);
      }
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
          // A cheap residue/tag check can avoid theoretical-spectrum generation.
          // Existing modification placement and cleavage-dependent digestion are
          // deliberately shared with enumeration in this first graph stage.
          const bool prefilter_rejected = graph_filter && std::none_of(low_it, up_it, [&](const auto& precursor) {
            return RNATagCandidateFilter::isCompatible(candidate, spectrum_graphs[precursor.second.scan_index], selected_ions,
                                                       precursor.second.charge);
          });
          if (prefilter_rejected)
          {
            ++graph_prefiltered;
            // Verification still builds the theory to check the rejection below.
            if (! verify_graph_matches) { continue; }
          }
          OPENMS_LOG_DEBUG << "Candidate: " << candidate.toString() << " ("
                           << float(candidate_mass) << " Da)" << endl;

          // pre-generate spectra:
          map<Int, MSSpectrum> theo_spectra_by_charge;
          spectrum_generator.getMultipleSpectra(theo_spectra_by_charge,
                                                candidate, precursor_charges,
                                                base_charge);

          if (!blacklist_mz.empty())
          {
            for (auto& [charge, spectrum] : theo_spectra_by_charge)
            {
              filterBlacklistedIons_(spectrum, blacklist_mz, blacklist_tolerance);
            }
          }
          map<Int, RNATagCandidateFilter::PreparedQuery> graph_queries;
          if (graph_enabled || (scoring_method == "route"))
          {
            for (const auto& [charge, model] : theo_spectra_by_charge)
            {
              graph_queries.emplace(charge, RNATagCandidateFilter::prepare(candidate, graph_alphabet, model, selected_ions, route_base_loss));
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
            if (diagnostic_filter && ! diagnostic_evidence[scan_index].accepts(candidate)) { continue; }
            if (prefilter_rejected)
            {
              // Verification only: the prefilter must not drop a pair that the
              // reference matcher supports or keeps through enumeration fallback.
              ++graph_verified;
              graph_mismatches += RNATagCandidateFilter::matchReference(candidate, spectrum_graphs[scan_index], exp_spectrum, theo_spectrum,
                                                                        selected_ions, search_param.fragment_mass_tolerance,
                                                                        search_param.fragment_tolerance_ppm)
                                    .accepted();
              continue;
            }
            RNATagCandidateFilter::Evidence graph_evidence;
            if (graph_enabled)
            {
              graph_evidence = RNATagCandidateFilter::match(graph_queries.at(charge * base_charge), spectrum_graphs[scan_index], exp_spectrum,
                                                            search_param.fragment_mass_tolerance, search_param.fragment_tolerance_ppm);
              if (verify_graph_matches)
              {
                const auto reference
                  = RNATagCandidateFilter::matchReference(candidate, spectrum_graphs[scan_index], exp_spectrum, theo_spectrum, selected_ions,
                                                          search_param.fragment_mass_tolerance, search_param.fragment_tolerance_ppm);
                ++graph_verified;
                graph_mismatches += graph_evidence != reference;
                graph_evidence = reference;
              }
              if (graph_evidence.supported)
              {
                ++graph_supported;
              }
              else if (! graph_evidence.fallback_reason.empty())
              {
                ++graph_fallback;
              }
              else
              {
                ++graph_rejected;
                if (graph_filter)
                {
                  continue;
                }
              }
            }
            vector<PeptideHit::PeakAnnotation> annotations;
            // Compute score using selected method
            double score;
            RNARouteScorer::Result route_evidence;
            if (scoring_method == "route")
            {
              if (route_score_limit_exceeded.load(std::memory_order_relaxed)) { continue; }
              try
              {
                route_evidence = RNARouteScorer::score(graph_queries.at(charge * base_charge), route_scoring_spectra[scan_index],
                                                       route_scoring_contexts[scan_index], negative_mode, route_score_options);
              }
              catch (const Exception::BaseException& error)
              {
                // Name the pair that failed: the limit is per model query, and
                // any other scoring error must not be reported as exhaustion.
#pragma omp critical (route_score_failure)
                {
                  if (route_score_failure.empty())
                  {
                    route_score_failure = "spectrum '" + exp_spectrum.getNativeID() + "', candidate '" + candidate.toString() + "': " + error.what();
                  }
                }
                route_score_limit_exceeded.store(true, std::memory_order_relaxed);
                continue;
              }
              score = route_evidence.score;
              annotations = route_evidence.annotations;
            }
            else
              {
                score = scoreNativeSpectrum_(scoring_method, search_param.fragment_mass_tolerance,
                search_param.fragment_tolerance_ppm, exp_spectrum, theo_spectrum,
                annotations, num_intensity_classes, tic_fraction, use_mass_accuracy, single_peak_credit);
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

            ++hit_counter;

            OPENMS_LOG_DEBUG << "Score: " << score << endl;

#pragma omp critical (annotated_hits_access)
            {
              HitsByScore& scan_hits = annotated_hits[scan_index];
              if (search_top_hits == 0 || scan_hits.size() < search_top_hits || score >= scan_hits.rbegin()->first)
              {
                AnnotatedHit new_hit;
                new_hit.oligo_ref = oligo_ref;
                new_hit.sequence = candidate;
                new_hit.precursor_error_ppm =
                  (prec_it->first - candidate_mass) / candidate_mass * 1.0e6;
                new_hit.annotations = annotations;
                new_hit.precursor_ref = &(prec_it->second);
                new_hit.graph_evaluated = graph_enabled;
                new_hit.graph_evidence = graph_evidence;
                new_hit.route_evidence = std::move(route_evidence);
                retainTopHit_(scan_hits, score, std::move(new_hit), search_top_hits);
              }
            }
          }
        }
      }
    }
    }
    progresslogger.endProgress();

    if (route_score_limit_exceeded.load(std::memory_order_relaxed))
    {
      OPENMS_LOG_ERROR << "RNA route scoring failed for " << route_score_failure << ". No search results were exported.\n";
      return UNEXPECTED_RESULT;
    }

    if (verify_graph_matches)
    {
      OPENMS_LOG_INFO << "RNA graph verified pairs: " << graph_verified << "; mismatches: " << graph_mismatches << '\n';
      if (graph_mismatches)
      {
        OPENMS_LOG_ERROR << "Graph prefiltering or prepared matching disagrees with the reference matcher.\n";
        return UNEXPECTED_RESULT;
      }
    }
    if (graph_generate)
    {
      OPENMS_LOG_INFO << "RNA placement nodes: " << placement_nodes << "; pruned branches: " << placement_pruned
                      << "; emitted variants: " << placement_emitted << "; enumeration fallback oligos: " << placement_fallbacks << '\n';
      if (placement_audit)
      {
        OPENMS_LOG_INFO << "RNA placement audit reference variants: " << placement_reference << "; eligible variants: " << placement_expected
                        << "; missing variants: " << placement_missing << "; extra variants: " << placement_extra << '\n';
        if (placement_missing || placement_extra)
        {
          OPENMS_LOG_ERROR << "Constrained modification placement disagrees with exhaustive graph filtering.\n";
          return UNEXPECTED_RESULT;
        }
      }
    }
    if (graph_enabled)
    {
      OPENMS_LOG_INFO << "RNA graph prefiltered variants: " << graph_prefiltered << "; supported pairs: " << graph_supported
                      << "; enumeration fallback pairs: " << graph_fallback << "; rejected pairs after theoretical validation: " << graph_rejected
                      << '\n';
    }
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
    const bool rescored =
      postProcessHits_(spectra, annotated_hits, id_data, negative_mode, file_ref, search_score_ref, report_top_hits, rescoring_scores);
    progresslogger.endProgress();
    OPENMS_LOG_INFO << "Identified spectra: " << id_data.getObservations().size()
                    << endl;

    // FDR:
    if (!decoy_pattern.empty())
    {
      OPENMS_LOG_INFO << "Performing FDR calculations..." << endl;
      calculateAndFilterFDR_(id_data, report_top_hits == 1, rescored ? rescoring_scores->score_ref : search_score_ref);
    }
    id_data.calculateCoverages();

    // store results
    if (!db_out.empty())
    {
      OMSFile(log_type_).store(db_out, id_data);
    }

    MzTab results = IdentificationDataConverter::exportMzTab(id_data,
      getFlag_("report:precursor_percentage") ? set<std::string>{"precursor_percentage"} : set<std::string>{},
      {"isotope_offset"});
    completeRNAScoreColumns_(results);
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
      if (!in_digest.empty())
      {
        // RNA seqs. were imported from a previous search run - need to "tag"
        // them with the current processing step so they get exported properly:
        for (IdentificationData::ParentSequenceRef ref =
               id_data.getParentSequences().begin(); ref !=
               id_data.getParentSequences().end(); ++ref)
        {
          // @TODO: find a way to avoid the copying (modify in place?):
          IdentificationData::ParentSequence copy = *ref;
          copy.addProcessingStep(id_data.getCurrentProcessingStep());
          id_data.registerParentSequence(copy);
        }
      }
      vector<ProteinIdentification> proteins;
      PeptideIdentificationList peptides;
      IdentificationDataConverter::exportIDs(id_data, proteins, peptides);
      FileHandler().storeIdentifications(id_out, proteins, peptides, {FileTypes::IDXML});
    }

    if (!lfq_out.empty())
    {
      generateLFQInput_(id_data, lfq_out);
    }

    if (!bedrmod_out.empty())
    {
      map<std::string, std::string> reference_metadata;
      for (const std::string key : {"organism", "assembly", "annotation_source", "annotation_version"})
      {
        reference_metadata[key] = getStringOption_("bedrmod:" + key);
      }
      reference_metadata["sequencing_platform"] = "ddMS2";
      reference_metadata["basecalling"] = "NASE";
      reference_metadata["bioinformatics_workflow"] = "NASE; caller-selected FDR cutoff and report:top_hits; score=min q-value (lower is better)";
      BedRModFile().store(bedrmod_out, id_data, bedrmod_chebi_mapping, reference_metadata);
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
