// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>

namespace OpenMS
{

/**
  @brief Minimal in-memory peptide-spectrum search engine.

  Searches MS2 spectra against a protein FASTA database and produces protein- and
  peptide-level identifications. Designed as a self-contained reference / teaching
  implementation; it is not intended as a feature-complete replacement for external
  engines such as MSGF+, Comet, Sage, or MSFragger.

  The pipeline run by search():
    -# Load spectra from @p in_spectra and the protein database from @p in_db.
    -# Preprocess MS2 spectra (filtering, deisotoping, decharging) via
       preprocessSpectra_() using the configured fragment tolerance.
    -# In-silico digest each protein with the configured enzyme/specificity/missed
       cleavages, generate all combinatorial fixed/variable modification variants
       (capped at @c modifications:variable_max_per_peptide), and score them
       against MS2 spectra whose precursor m/z matches the candidate peptide mass
       within @c precursor:mass_tolerance (with optional 1/-1 isotope correction
       per @c precursor:isotopes).
    -# Collect the top-N scoring candidates per spectrum (@c report:top_hits),
       optionally generate target/decoy results (@c decoys), filter by FDR
       (@c FDR:PSM, q-value) and annotate PSMs (@c annotate:PSM) in
       postProcessHits_().

  Configuration is exposed through the @ref DefaultParamHandler base — see the
  defaults installed by the constructor for the full list of supported keys
  (precursor/fragment tolerances and units, charge range, modifications, enzyme,
  peptide size/motif/missed-cleavage filters, FDR threshold, top-hits, etc.).

  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI SimpleSearchEngineAlgorithm :
  public DefaultParamHandler,
  public ProgressLogger
{
  public:
    /// Default constructor; installs the search parameters (see class docs)
    SimpleSearchEngineAlgorithm();

    /// Outcome of search(), distinguishing recoverable input issues from execution errors
    enum class ExitCodes
    {
      EXECUTION_OK,       ///< Search completed; @p prot_ids and @p pep_ids contain the result
      INPUT_FILE_EMPTY,   ///< Spectrum input contained no usable MS2 spectra after loading/filtering
      UNEXPECTED_RESULT,  ///< Internal post-condition violated (e.g. no candidates scored at all)
      UNKNOWN_ERROR,      ///< Caught a generic exception; details written to the log
      ILLEGAL_PARAMETERS  ///< Configuration is internally inconsistent or unsupported
    };

    /**
      @brief Search the MS2 spectra in @p in_spectra against the protein database in @p in_db.

      Spectra and database are loaded from disk; the result is written into the two
      output arguments. Existing contents of @p prot_ids and @p pep_ids are not
      cleared by this call. The current parameter set (see the class brief) controls
      tolerances, modifications, enzyme, FDR, etc.

      @param[in]  in_spectra Path to the spectrum input (mzML or any format readable by @ref FileHandler).
      @param[in]  in_db      Path to the protein FASTA database to search against.
      @param[out] prot_ids   Protein identifications produced by the search (one run per call).
      @param[out] pep_ids    Peptide-spectrum matches (PSMs) produced by the search.
      @return Status code; see @ref ExitCodes.
    */
    ExitCodes search(const String& in_spectra,
      const String& in_db,
      std::vector<ProteinIdentification>& prot_ids,
      PeptideIdentificationList& pep_ids) const;
  protected:
    void updateMembers_() override;

    /**
      @brief Compact internal record for one scored peptide candidate against one spectrum.

      Used instead of @ref PeptideHit while scoring so that all candidates per spectrum
      fit in memory; only the top-N (@c report:top_hits) entries per spectrum are later
      promoted into @ref PeptideHit objects by postProcessHits_(). Field order is chosen
      to minimise struct padding (doubles first, then floats, then ints).
    */
    struct AnnotatedHit_
    {
      StringView sequence; ///< unmodified peptide sequence (view into the digested protein)
      SignedSize peptide_mod_index; ///< enumeration index of the modification variant (for re-materialisation)
      // Layout: doubles first, then floats, then int, then uint16_t — minimizes padding
      double score = 0; ///< main score (higher is better)
      float prefix_fraction = 0; ///< fraction of annotated prefix ions (a/b/c)
      float suffix_fraction = 0; ///< fraction of annotated suffix ions (x/y/z)
      float mean_error = 0.0f; ///< mean absolute fragment mass error
      int isotope_error = 0; ///< matched precursor isotope offset (per @c precursor:isotopes)
      uint16_t matched_prefix_ions = 0; ///< number of matched prefix ions (a/b/c)
      uint16_t matched_suffix_ions = 0; ///< number of matched suffix ions (x/y/z)

      /// Strict-weak ordering: higher score first, then lower mod index, then sequence
      static bool hasBetterScore(const AnnotatedHit_& a, const AnnotatedHit_& b)
      {
        if (a.score != b.score) return a.score > b.score;
        if (b.peptide_mod_index != a.peptide_mod_index) return a.peptide_mod_index < b.peptide_mod_index;
        return a.sequence < b.sequence;
      }
    };

    /**
      @brief Preprocess MS2 spectra in place: filter, deisotope, decharge.

      Applies the standard search-engine spectrum normalisation used before scoring:
      filtering out low-quality peaks, charge state deconvolution, and isotope-pattern
      deisotoping using the supplied fragment tolerance.

      @param[in,out] exp                              Spectra to preprocess in place.
      @param[in]     fragment_mass_tolerance          Tolerance for deisotoping and decharging.
      @param[in]     fragment_mass_tolerance_unit_ppm If true, @p fragment_mass_tolerance is ppm; otherwise Th.
    */
    static void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm);

    /**
      @brief Materialise top-N scored candidates per spectrum into @ref PeptideHit / @ref ProteinIdentification objects.

      Converts the in-memory @ref AnnotatedHit_ records produced by the scoring loop into
      first-class identification objects: re-applies the modification variant indicated
      by @c peptide_mod_index, annotates PSM meta values (per @c annotate:PSM), populates
      protein references, and stamps search-engine settings (tolerances, modifications,
      enzyme, etc.) onto the resulting @ref ProteinIdentification so the output is
      self-describing.

      Most parameters mirror the algorithm's own configuration and are passed in
      explicitly so this routine can also be reused outside member context.

      @param[in]  exp                              Preprocessed spectra used as scoring input.
      @param[in]  annotated_hits                   Per-spectrum vectors of scored candidates (already top-N filtered).
      @param[out] protein_ids                      Protein identifications to populate.
      @param[out] peptide_ids                      Peptide-spectrum matches to populate.
      @param[in]  top_hits                         Maximum number of PSMs per spectrum to materialise.
      @param[in]  fixed_modifications              Resolved fixed-modification table.
      @param[in]  variable_modifications           Resolved variable-modification table.
      @param[in]  max_variable_mods_per_peptide    Cap on simultaneous variable modifications.
      @param[in]  modifications_fixed              UniMod names of fixed modifications (for the ID metadata stamp).
      @param[in]  modifications_variable           UniMod names of variable modifications (for the ID metadata stamp).
      @param[in]  peptide_missed_cleavages         Allowed missed cleavages (for the ID metadata stamp).
      @param[in]  precursor_mass_tolerance         Precursor mass tolerance value.
      @param[in]  fragment_mass_tolerance          Fragment mass tolerance value.
      @param[in]  precursor_mass_tolerance_unit_ppm "ppm" or "Da"; recorded in the ID metadata.
      @param[in]  fragment_mass_tolerance_unit_ppm  "ppm" or "Da"; recorded in the ID metadata.
      @param[in]  precursor_min_charge             Minimum precursor charge considered.
      @param[in]  precursor_max_charge             Maximum precursor charge considered.
      @param[in]  enzyme                           Enzyme name (for the ID metadata stamp).
      @param[in]  database_name                    FASTA database path/name (for the ID metadata stamp).
    */
    void postProcessHits_(const PeakMap& exp,
      std::vector<std::vector<SimpleSearchEngineAlgorithm::AnnotatedHit_> >& annotated_hits,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      Size top_hits,
      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
      Size max_variable_mods_per_peptide,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const String& precursor_mass_tolerance_unit_ppm,
      const String& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const String& enzyme,
      const String& database_name) const;

    double precursor_mass_tolerance_;       ///< Precursor mass tolerance (value); unit in @c precursor_mass_tolerance_unit_
    String precursor_mass_tolerance_unit_;  ///< "ppm" or "Da"

    Size precursor_min_charge_;             ///< Minimum precursor charge considered
    Size precursor_max_charge_;             ///< Maximum precursor charge considered

    IntList precursor_isotopes_;            ///< Allowed precursor isotope offsets (0 = monoisotopic, 1 = +1 Da, etc.)

    double fragment_mass_tolerance_;        ///< Fragment mass tolerance (value); unit in @c fragment_mass_tolerance_unit_

    String fragment_mass_tolerance_unit_;   ///< "ppm" or "Da"

    StringList modifications_fixed_;        ///< UniMod names of fixed modifications

    StringList modifications_variable_;     ///< UniMod names of variable modifications

    Size modifications_max_variable_mods_per_peptide_; ///< Cap on simultaneous variable modifications per peptide

    String enzyme_;                         ///< Enzyme name as recognised by @ref EnzymaticDigestion

    bool decoys_;                           ///< If true, generate target/decoy results

    double fdr_psm_;                        ///< q-value threshold for PSM filtering (0 = disabled); requires @c decoys_

    StringList annotate_psm_;               ///< PSM meta-value annotations to add (see @c annotate:PSM defaults)

    Size peptide_min_size_;                 ///< Minimum peptide length after digestion
    Size peptide_max_size_;                 ///< Maximum peptide length after digestion (0 = unlimited)
    Size peptide_missed_cleavages_;         ///< Allowed missed cleavages in digestion
    EnzymaticDigestion::Specificity peptide_enzyme_specificity_{EnzymaticDigestion::SPEC_FULL}; ///< full / semi / none

    String peptide_motif_;                  ///< Optional regex motif; only peptides matching are considered

    Size report_top_hits_;                  ///< Number of top-scoring PSMs reported per spectrum
};

} // namespace

