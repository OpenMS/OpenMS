// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>


#include <array>
#include <mutex>
#include <vector>
#include <functional>
#include <algorithm>   // std::max (used by inline static isOpenSearchMode)

namespace OpenMS
{
  /** @brief Generates from a set of Fasta files a 2D-datastructure which stores all theoretical masses of all
   * b and y ions from all peptides generated from the Fasta file. The datastructure is build such that on one axis
   * the fragments are sorted by their own mass and the axis by the mass of their precursor/protein.
   * The FI has two options: Bottom-up and Top Down. In later digestion is skiped and the fragments have a direct
   * reference to the mass of the proteins instead of digested peptides.
   */
  class OPENMS_DLLAPI FragmentIndex : public DefaultParamHandler
  {
  public:


    /** @brief Compact descriptor of a peptide instance held by the FragmentIndex.
     *
     * Field semantics and how they relate to the braced initializer lists used in tests {a, b, {c, d}, e}:
     *  - protein_idx .......... 'a' Index into the FASTA entries vector passed to build(); identifies the source protein.
     *  - mod_bitmask_ ......... 'b' Bitmask encoding which variable modification slots are active.
     *                             Each bit corresponds to a (position, mod_type) pair found by scanning
     *                             the sequence left-to-right against the variable modification config.
     *                             0 = unmodified (or fixed-only). Supports up to 32 variable mod slots.
     *  - sequence_ ............ '{c, d}' 0-based start offset and length (in residues) of the peptide subsequence within the protein.
     *                             Note: length (d) is used like std::string::substr(start, length).
     *  - precursor_mz_ ........ 'e' Mono-isotopic m/z at charge 1 (M+H)+; used for ordering/slicing in the index.
     *                             Many tests use a dummy value here since only ordering invariants are asserted.
     */
    struct Peptide {

      // We need a constructor in order to emplace back
      Peptide(UInt32 protein_idx, uint32_t mod_bitmask, std::pair<uint16_t , uint16_t> sequence, float precursor_mz):
          protein_idx(protein_idx),
        mod_bitmask_(mod_bitmask),
        sequence_(sequence),
        precursor_mz_(precursor_mz)
        {}

        UInt32 protein_idx;            ///< 0-based index into FASTA entries provided to build(); identifies the source protein
        uint32_t mod_bitmask_;         ///< Bitmask of active variable mod slots (0 = unmodified/fixed-only; up to 32 slots)
        std::pair<uint16_t , uint16_t> sequence_; ///< {start, length} within the source protein sequence (start is 0-based; length in residues)
        float precursor_mz_;           ///< Mono-isotopic m/z at charge 1 (M+H)+ of this peptide; used for sorting/filtering
    };

    /**
     * @brief Match between a query peak and an entry in the DB
     */
    struct SpectrumMatch
    {
      uint32_t num_matched_{};       ///< Number of peaks-fragment hits
      uint32_t subset_bitmask_{};    ///< SNES v1.1: active slots in the slot list returned by buildModSlots_. 0 = unmodified. Ignored in non-SNES mode.
      float    sigma_delta_{};       ///< SNES v1.1: Σ of variable-mod deltas for this match. 0 in non-SNES / unmodified SNES.
      uint16_t precursor_charge_{};  ///< The precursor_charge used for the performed search
      int16_t  isotope_error_{};     ///< The isotope_error used for the performed search
      size_t   peptide_idx_{};       ///< The idx this struct belongs to
    };


    /**
     * @brief container for SpectrumMatch. Also keeps count of total number of candidates and total number of matches.
     */
    struct SpectrumMatchesTopN
    {
      std::vector<SpectrumMatch> hits_;     ///< The preliminary candidates


      SpectrumMatchesTopN() = default;

      /**
       * @brief Appends the a SpectrumMatchesTopN to another one. Add the number of all matched peaks up. Same for number of scored candidates
       * The
       * @param[in] other The appended struct
       * @return The struct after the attachment
       */
      SpectrumMatchesTopN& operator+=(const SpectrumMatchesTopN& other)
      {

        this->hits_.insert(this->hits_.end(), other.hits_.begin(), other.hits_.end());
        return *this;
      }

      void clear()
      {
        hits_.clear();

      }
    };
    /**
     * @brief Default constructor.
     *
     * Initializes an empty FragmentIndex. Call build() before using any query
     * functions. After clear(), the index returns to this unbuilt state.
     *
     * Thread-safety: constructing the object is thread-safe as long as the instance
     * is not shared across threads before initialization completes.
     */
    FragmentIndex();

    /**
     * @brief Default destructor.
     *
     * Releases owned memory. If the index was built, all internal buffers and
     * fragment buckets are freed. No exceptions are thrown.
     */
    ~FragmentIndex() override = default;

    /**
     * @brief Indicates whether the fragment index has been built.
     *
     * @return true if build() has completed successfully and the index is ready
     *         for queries; false otherwise (e.g., after construction or after clear()).
     *
     * Thread-safety: read-only and can be called concurrently with other
     * read-only methods. Must not race with build()/clear() on the same instance.
     */
    bool isBuild() const;

    /**
     * @brief Returns a reference to the internal peptide container.
     *
     * Provides read-only access to all peptides currently held by the index,
     * typically populated during build().
     *
     * @return const reference to the internal std::vector of Peptide.
     *
     * Preconditions: The vector may be empty if build() has not been called yet.
     * Thread-safety: read-only view; safe to access concurrently as long as no
     * thread mutates the index (e.g., build()/clear()).
     */
    const std::vector<Peptide>& getPeptides() const;

#ifdef DEBUG_FRAGMENT_INDEX
    /**
     * @brief Manually adds a peptide to the internal peptide list (debug builds only).
     *
     * Allows injecting a custom peptide sequence into the index prior to building,
     * e.g., for targeted testing. This function modifies the internal state and
     * must be used with care.
     *
     * @param[in,out] peptide AASequence of the peptide to add. The sequence may be modified
     *                internally (e.g., normalization/annotation steps).
     * @param[in] source_idx Index of the originating FASTA entry (or synthetic source)
     *                   to maintain provenance in downstream processing.
     *
     * Preconditions:
     *  - Must be called after peptides have been generated (e.g., generatePeptides())
     *    and before build(). Calling it after build() leads to undefined behavior.
     *
     * Thread-safety:
     *  - Not thread-safe. Do not call concurrently with build(), clear(), or any
     *    read operations. Restrict usage to single-threaded setup in debug builds.
     *
     * Exceptions:
     *  - Strong exception guarantee: either the peptide is added or the index remains unchanged.
     */
    void addSpecialPeptide(AASequence& peptide, Size source_idx);
#endif

    /** @brief Given a set of Fasta files, builds the Fragment Index datastructure (FID). First all fragments are sorted
     * by their own mass. Next they are placed in buckets. The min-fragment mass is stored for each bucket, whereupon
     * the fragments are sorted within the buckets by their originating precursor mass.
     *
     * @param[in] fasta_entries The FASTA entries used to build the index.
     */
    void build(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);

    /** @brief Delete fragment index. Sets is_build=false*/
    void clear();


    /** Return the [begin_idx, end_idx) peptide index range such that
     * `fi_peptides_[i].precursor_mz_ ∈ [precursor_mass + window.first, precursor_mass + window.second]`
     * for all i in the returned range.
     *
     * @param[in] precursor_mass The mono-charged precursor mass (M+H).
     * @param[in] window Signed absolute offsets around the precursor mass. By convention
     *                   `window.first` is <= 0 and `window.second` is >= 0 (produced by
     *                   `computeMassWindow_`). A reversed window trivially returns an empty
     *                   range; no diagnostic is emitted. No hidden tolerance is added.
     * @return [begin_idx, end_idx) half-open index range into `fi_peptides_`.
     */
    std::pair<size_t, size_t> getPeptidesInMassWindow(float precursor_mass,
                                                      const std::pair<float, float>& window) const;

    /// Shared auto-detection: open-search iff max(lower, upper) > threshold (1000 ppm or 1 Da).
    /// Strict `>`: exactly 1000 ppm stays closed.
    /// This is the single source of truth for the open-search auto-detection rule and is
    /// reused by ProSEAlgorithm and the TOPP tool.
    static bool isOpenSearchMode(double lower_magnitude,
                                 double upper_magnitude,
                                 bool unit_ppm) noexcept
    {
      const double threshold = unit_ppm ? 1000.0 : 1.0;
      return std::max(lower_magnitude, upper_magnitude) > threshold;
    }

    /// @name SNES (Speedy Non-specific Enzyme Search) bit encoding
    ///
    /// When the index is built in SNES mode (@ref isSnesMode), a @ref Peptide entry
    /// represents a *mother peptide* — the longest peptide anchored at one terminus of
    /// the protein. Bit 31 of @c Peptide::mod_bitmask_ distinguishes Single-N mothers
    /// (N-terminus anchored, C-terminus free, b-ion series indexed) from Single-C
    /// mothers (C-terminus anchored, N-terminus free, y-ion series indexed). SNES v1
    /// does not enumerate variable modifications on mothers, so the lower 31 bits are
    /// currently unused (reserved for variable-mod-on-mother support in a future
    /// version).
    ///
    /// When the index is built in the default (non-SNES) mode, @c mod_bitmask_ uses
    /// all 32 bits as variable-modification slots (existing semantics) and these
    /// accessors are not consulted.
    ///
    /// Rationale for stealing a single bit from the bitmask instead of adding a field
    /// to @c Peptide: zero runtime overhead and zero size impact for the common
    /// full-specificity trypsin case — the bit is only interpreted when the index
    /// dispatcher has already branched on @c is_snes_mode_.
    static constexpr uint32_t SNES_KIND_BIT_MASK = 1u << 31;   ///< bit 31; set = Single-C mother
    static constexpr uint32_t SNES_SLOT_MASK = ~SNES_KIND_BIT_MASK; ///< bits 0..30 in SNES mode

    /// SNES v1.1: constrain bin-walk hits to mothers with a specific protein
    /// anchor. Used to gate walks that enumerate PROTEIN_N_TERM / PROTEIN_C_TERM
    /// variable mods.
    enum class SnesAnchor
    {
      NONE,        ///< no anchor restriction (baseline walks)
      PROT_NTERM,  ///< mother must have sequence_.first == 0
      PROT_CTERM   ///< mother must have sequence_.first + sequence_.second == protein length
    };

    /// @return true iff the mother described by @p mod_bitmask is a Single-C (C-anchored) mother.
    /// Only meaningful for peptides from an SNES-built index.
    static bool isSingleCMother(uint32_t mod_bitmask) noexcept
    {
      return (mod_bitmask & SNES_KIND_BIT_MASK) != 0;
    }
    /// @return true iff the mother described by @p mod_bitmask is a Single-N (N-anchored) mother.
    static bool isSingleNMother(uint32_t mod_bitmask) noexcept
    {
      return (mod_bitmask & SNES_KIND_BIT_MASK) == 0;
    }

    /// @return true if the index was built in SNES mode.
    /// SNES activates when @em both @c snes_enabled is set to true and
    /// @c peptide:enzyme_specificity is @c none. @c snes_enabled defaults to false
    /// in v1 (opt-in), so specific/semi-specific searches and non-specific searches
    /// without the flag produce the same fragment index as the pre-SNES code path.
    bool isSnesMode() const noexcept { return is_snes_mode_; }


    /**
     * A match between a single query peak and a database fragment
     */
    struct Hit
    {
      Hit(UInt32 peptide_idx, float fragment_mz) :
        peptide_idx(peptide_idx),
        fragment_mz(fragment_mz)
      {}
      UInt32 peptide_idx; // index in database
      float fragment_mz;
    };

    /**@brief Queries one peak
     * @param[in] peak The queried peak
     * @param[in] peptide_idx_range The range of precursors/peptides the peptide could potentially belongs to
     * @param[in] peak_charge The charge of the peak. Is used to calculate the mass from the mz
     * @return a vector of Hits(matching peptide_idx_range and matching fragment_mz_) containing the idx of the hitted peptide and the mass of the hit
     */
    std::vector<Hit> query(const Peak1D& peak,
                           const std::pair<size_t,size_t>& peptide_idx_range,
                           uint16_t peak_charge);

    /**
     * @brief: queries one complete experimental spectra against the Database. Loops over all precursor charges
     * Starts at min_precursor_charge and iteratively goes to max_precursor_charge. We query all peaks multiple times with all the
     * different precursor charges and corresponding precursor masses
     * @param[in] spectrum experimental spectrum
     * @param[out] sms The n best Spectrum matches
     */
    void querySpectrum(const MSSpectrum& spectrum,
                       SpectrumMatchesTopN& sms);

    /**
     * @brief Query a spectrum against the fragment index with FASTA context.
     *
     * Required when FragmentIndex is in SNES mode with variable modifications;
     * the FASTA is needed to realize sub-peptide sequences and apply variable mods.
     * Non-SNES and SNES-without-var-mods paths ignore the @p fasta_entries argument.
     *
     * @param[in]  spectrum      Experimental spectrum with a single precursor.
     * @param[in]  fasta_entries The FASTA database passed to build().
     * @param[out] sms           Accumulated candidate matches.
     */
    void querySpectrum(const MSSpectrum& spectrum,
                       const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                       SpectrumMatchesTopN& sms);

    /** @brief Reconstruct a fully modified AASequence from a Peptide's bitmask.
     *
     * Used for result output - only called for final hits (not in the build hot path).
     * Applies fixed modifications, then uses the bitmask to determine which variable
     * modifications are active at which positions.
     *
     * @param[in] peptide   The Peptide descriptor with mod_bitmask_
     * @param[in] fasta_entries The FASTA database used during build()
     * @return The fully modified AASequence
     */
    AASequence reconstructModifiedSequence(const Peptide& peptide,
                                           const std::vector<FASTAFile::FASTAEntry>& fasta_entries) const;

    /** @brief Find the realized sub-peptide length of a SNES mother that best matches
     * the observed precursor mass.
     *
     * Scans realizable lengths k in [peptide_min_length_, mother.sequence_.second] from
     * the appropriate terminus (left-to-right for Single-N mothers, right-to-left for
     * Single-C mothers), computing the cumulative residue mass plus fixed modifications.
     * Returns the length whose realized (M+H)+ mass is closest to @p target_mh_plus
     * within the given tolerance, or -1 if no length satisfies the tolerance.
     *
     * Must only be called in SNES mode; returns -1 immediately otherwise.
     *
     * @param[in] mother A Peptide representing a Single-N or Single-C mother.
     * @param[in] fasta_entries Source database (same one passed to build()).
     * @param[in] target_mh_plus Observed (M+H)+ mass after isotope-error correction.
     * @param[in] tolerance_lower_magnitude Positive tolerance magnitude on the low side
     *            (realized_mass - target >= -tolerance_lower_magnitude).
     * @param[in] tolerance_upper_magnitude Positive tolerance magnitude on the high side
     *            (realized_mass - target <= +tolerance_upper_magnitude).
     * @param[in] tolerance_ppm True if the magnitudes are in ppm.
     * @return Realized length, or -1 on failure.
     */
    int realizeSNESLength(const Peptide& mother,
                          const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                          double target_mh_plus,
                          double tolerance_lower_magnitude,
                          double tolerance_upper_magnitude,
                          bool tolerance_ppm) const;

    /// Reconstruct a realized SNES sub-peptide as an AASequence.
    /// @param mother the SNES mother Peptide entry
    /// @param fasta_entries the FASTA entries used to build the index
    /// @param realized_length the length of the realized sub-peptide (from realizeSNESLength)
    /// @param subset_bitmask SNES v1.1: active slots from buildModSlots_(seq_ptr, realized_length, ...)
    ///        to apply as variable modifications. 0 = unmodified (backward compatible).
    /// @return AASequence representing the realized sub-peptide with any variable mods from subset_bitmask applied
    AASequence reconstructRealizedSubSequence(const Peptide& mother,
                                              const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                              size_t realized_length,
                                              uint32_t subset_bitmask = 0) const;

protected:


  /**@brief One entry in the fragment index
   */
  struct Fragment
  {
      Fragment() = default;
      Fragment(UInt32 peptide_idx, float fragment_mz):
          peptide_idx_(peptide_idx),
          fragment_mz_(fragment_mz)
      {}
      UInt32 peptide_idx_{}; // 32 bit in sage
      float fragment_mz_{};
  };

    bool is_build_{false};              ///< true, if the database has been populated with fragments

    void updateMembers_() override;

     /**@brief Generates all peptides from given fasta entries. If Bottom-up is set to false
     * skips digestion. If set to true the Digestion enzyme can be set in the parameters.
     * Additionally introduces fixed and variable modifications for restrictive PSM search.
     *
     * @param[in] fasta_entries
     */
    void generatePeptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

    /**@brief SNES-mode peptide enumeration: emit Single-N + Single-C mother peptides.
     *
     * Called instead of the usual enzymatic-digestion path when @c is_snes_mode_ is
     * true. For each protein, emits:
     *  - one Single-N mother at every anchor position i in [0, L - min_length], whose
     *    sequence spans residues [i, i + min(max_length, L - i)). Tagged by clearing
     *    bit @c SNES_KIND_BIT_MASK in the mod_bitmask; indexed with b-ion series only.
     *  - one Single-C mother at every anchor position j in [min_length - 1, L - 1],
     *    whose sequence spans residues [j - min(max_length, j + 1) + 1, j + 1). Tagged
     *    by setting bit @c SNES_KIND_BIT_MASK; indexed with y-ion series only.
     *
     * All sub-peptides of the mother (down to the configured min_length) can be
     * realized from a single mother record, which is what gives SNES its memory and
     * speed win over naïve O(L^2) non-specific enumeration.
     *
     * v1 restriction: variable modifications are disabled in SNES mode. A warning is
     * emitted at build time if any variable modification is configured. Fixed
     * modifications (both residue-specific and terminal) are fully supported.
     *
     * @param[in] fasta_entries  Protein database (same semantics as generatePeptides).
     */
    void generateSNESMothers_(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

    /** @brief Entry in the per-AA variable modification lookup table. */
    struct VarModEntry
    {
      double delta_mass;                    ///< mass delta from this modification
      const ResidueModification* mod_ptr;   ///< pointer to the modification (for AASequence reconstruction)
      ResidueModification::TermSpecificity term_spec; ///< where this mod can be applied
    };

    /** @brief A candidate modification slot for a specific peptide.
     *
     * Slots are built by scanning a peptide sequence left-to-right against the variable mod config.
     * The slot index determines its bit position in mod_bitmask_.
     */
    struct ModSlot
    {
      uint16_t position;                    ///< residue index, or NTERM_SLOT/CTERM_SLOT
      double delta_mass;                    ///< mass delta
      const ResidueModification* mod_ptr;   ///< for AASequence reconstruction

      static constexpr uint16_t NTERM_SLOT = UINT16_MAX - 1; ///< sentinel for pure N-terminal mod slot
      static constexpr uint16_t CTERM_SLOT = UINT16_MAX;      ///< sentinel for pure C-terminal mod slot
    };

    static constexpr size_t MAX_MOD_SLOTS = 32; ///< max variable mod slots per peptide (uint32_t bitmask)

    /// Build per-AA modification lookup tables from modifications_fixed_ and modifications_variable_.
    /// Called once at the start of generatePeptides().
    void initModificationTables_();

    /// Scan a peptide sequence to find all variable modification slots.
    /// Returns the number of slots written to out_slots (at most MAX_MOD_SLOTS).
    /// Deterministic ordering: N-term pure-terminal mods, then left-to-right residue mods
    /// (ANYWHERE + position-specific terminal), then C-term pure-terminal mods.
    /// @param sequence raw amino acid character array
    /// @param seq_len length of the sequence
    /// @param out_slots output array for modification slots (must have space for MAX_MOD_SLOTS entries)
    /// @param is_protein_nterm true if this peptide starts at protein position 0
    /// @param is_protein_cterm true if this peptide ends at the last protein residue
    size_t buildModSlots_(const char* sequence, size_t seq_len, ModSlot* out_slots,
                          bool is_protein_nterm = false, bool is_protein_cterm = false) const;

    /// Enumerate distinct Σ values achievable by any subset of configured
    /// variable mods with popcount ≤ max_variable_mods_per_peptide_.
    /// Configuration-global; does not consider per-peptide residue inventory.
    /// Per-peptide applicability is enforced at query-time subset enumeration.
    ///
    /// @param include_prot_nterm_mods include mods with PROTEIN_N_TERM specificity
    /// @param include_prot_cterm_mods include mods with PROTEIN_C_TERM specificity
    /// @return sorted ascending distinct Σ values; always includes 0.0
    std::vector<double> computeSnesSigmaDeltaSet_(bool include_prot_nterm_mods,
                                                   bool include_prot_cterm_mods) const;

    /// Per-AA fixed modification delta mass (0.0 if no fixed mod applies)
    std::array<double, 128> fixed_mod_deltas_{};
    /// Per-AA fixed modification pointer (nullptr if none)
    std::array<const ResidueModification*, 128> fixed_mod_ptrs_{};
    double fixed_nterm_delta_{0.0};   ///< Fixed N-terminal mod delta (0 if none)
    double fixed_cterm_delta_{0.0};   ///< Fixed C-terminal mod delta (0 if none)
    const ResidueModification* fixed_nterm_mod_ptr_{nullptr};
    const ResidueModification* fixed_cterm_mod_ptr_{nullptr};

    /// Per-AA variable modification table: for each ASCII char, list of possible variable mods
    std::array<std::vector<VarModEntry>, 128> variable_mod_table_{};
    /// Pure N-terminal variable mods (not residue-specific)
    std::vector<VarModEntry> variable_nterm_mods_;
    /// Pure C-terminal variable mods (not residue-specific)
    std::vector<VarModEntry> variable_cterm_mods_;

    bool mod_tables_initialized_{false};

    /// SNES mode state. Set in @c updateMembers_ from the @c snes_enabled parameter.
    /// When true, @ref generatePeptides dispatches to @ref generateSNESMothers_ and
    /// the fragment-index query layer switches to the one-sided lookup. When false,
    /// no SNES code path is active and the index behaves identically to the original
    /// precursor-window-based implementation.
    bool is_snes_mode_{false};

    /// User-facing SNES opt-in switch (parameter "snes_enabled"). Only takes effect
    /// when the configured enzyme specificity is @c SPEC_NONE — specific / semi-
    /// specific searches ignore it. Exposed as a separate member so the parameter
    /// can be set/queried independently of the derived @c is_snes_mode_ state
    /// (which captures the combined decision specificity && snes_enabled).
    bool snes_enabled_{false};

    /// SNES v1.1: precomputed distinct Σ_delta values for bin-walk targets.
    /// Baseline set excludes protein-term-only variable mods.
    std::vector<double> snes_sigma_delta_set_;
    /// SNES v1.1: Σ values including PROTEIN_N_TERM-only variable mods.
    /// Used only for Single-N mothers anchored at protein position 0.
    std::vector<double> snes_sigma_delta_set_with_prot_nterm_;
    /// SNES v1.1: Σ values including PROTEIN_C_TERM-only variable mods.
    /// Used only for Single-C mothers anchored at the protein C-terminus.
    std::vector<double> snes_sigma_delta_set_with_prot_cterm_;

    /// Precomputed residue mass lookup table: ASCII char -> internal monoisotopic mass (Da).
    /// Indexed by single-letter amino acid code (e.g., 'A'=65). Entries for non-AA chars are 0.
    static std::array<double, 128> residue_mass_table_;
    static std::once_flag mass_table_once_flag_;
    static void initResidueMassTable_();

    /// Precomputed ion-type mass offsets (from Residue::getInternalTo*Ion formulas)
    struct IonOffsets
    {
      double b_offset{0.0};
      double y_offset{0.0};
      double a_offset{0.0};
      double c_offset{0.0};
      double x_offset{0.0};
      double z_offset{0.0};
    };
    static IonOffsets ion_offsets_;

    /// Lightweight fragment generation: compute b/y ion m/z directly from amino acid chars.
    /// Bypasses AASequence::fromString and TheoreticalSpectrumGenerator.
    /// Uses the class-level @c add_b_ions_ / @c add_y_ions_ / ... flags for the ion
    /// series selection. See @ref generateFragmentsForSeries_ for the explicit-flag
    /// variant used by the SNES mother path.
    /// @param[out] fragments  Output vector to append Fragment entries to
    /// @param[in]  sequence   Raw amino acid string (no modifications)
    /// @param[in]  seq_len    Length of sequence
    /// @param[in]  peptide_idx Index of this peptide in fi_peptides_
    /// @param[in]  n_term_mod_mass  Mass delta from N-terminal modification (0 if none)
    /// @param[in]  c_term_mod_mass  Mass delta from C-terminal modification (0 if none)
    /// @param[in]  residue_mod_masses  Per-residue modification mass deltas (nullptr if none; array of seq_len doubles)
    void generateFragmentsLightweight_(
      std::vector<Fragment>& fragments,
      const char* sequence,
      size_t seq_len,
      UInt32 peptide_idx,
      double n_term_mod_mass,
      double c_term_mod_mass,
      const double* residue_mod_masses) const;

    /// Fragment generation with explicit per-call ion-series selection.
    ///
    /// Called by the SNES mother path to restrict a Single-N mother to b-ions and a
    /// Single-C mother to y-ions regardless of the class-level @c add_*_ions_ flags.
    /// @c generateFragmentsLightweight_ forwards to this function after packing the
    /// class flags; both share a single implementation.
    ///
    /// @param[out] fragments Output vector to append Fragment entries to
    /// @param[in] sequence Raw amino acid string (no modifications)
    /// @param[in] seq_len Length of sequence
    /// @param[in] peptide_idx Index of this peptide in fi_peptides_
    /// @param[in] n_term_mod_mass Mass delta from N-terminal modification (0 if none)
    /// @param[in] c_term_mod_mass Mass delta from C-terminal modification (0 if none)
    /// @param[in] residue_mod_masses Per-residue modification mass deltas (nullptr if none; array of seq_len doubles)
    /// @param[in] add_b Emit b-ions (prefix).
    /// @param[in] add_a Emit a-ions (prefix).
    /// @param[in] add_c Emit c-ions (prefix).
    /// @param[in] add_y Emit y-ions (suffix).
    /// @param[in] add_x Emit x-ions (suffix).
    /// @param[in] add_z Emit z-ions (suffix).
    void generateFragmentsForSeries_(
      std::vector<Fragment>& fragments,
      const char* sequence,
      size_t seq_len,
      UInt32 peptide_idx,
      double n_term_mod_mass,
      double c_term_mod_mass,
      const double* residue_mod_masses,
      bool add_b,
      bool add_a,
      bool add_c,
      bool add_y,
      bool add_x,
      bool add_z) const;

    std::vector<Peptide> fi_peptides_;   ///< vector of all (digested) peptides
    std::vector<Fragment> fi_fragments_; ///< vector of all theoretical fragments (b- and y- ions)

    /// Protein lengths indexed by protein_idx, populated at build() time.
    /// Used by SNES v1.1 to gate PROTEIN_C_TERM variable-mod bin walks.
    std::vector<uint32_t> protein_lengths_;

    float fragment_min_mz_;  ///< smallest fragment mz
    float fragment_max_mz_;  ///< largest fragment mz
    size_t min_ion_index_{0}; ///< skip ions below this index (0=all, 2=skip b1/b2/y1/y2)
    size_t bucketsize_;       ///< number of fragments per outer node
    std::vector<float> bucket_min_mz_;  ///< vector of the smalles fragment mz of each bucket
    double precursor_mass_tolerance_lower_{20.0};   ///< positive magnitude, effective lower bound is -lower
    double precursor_mass_tolerance_upper_{20.0};   ///< positive magnitude, effective upper bound is +upper
    bool precursor_mass_tolerance_unit_ppm_{true};
    float fragment_mz_tolerance_;
    bool fragment_mz_tolerance_unit_ppm_{true};    
private:


    /**
     * @brief SNES-mode spectrum query (MetaMorpheus-style: byte-count + b-ion filter).
     *
     * Implements the Rolfs/Smith 2020 inverted-index search strategy as executed in
     * MetaMorpheus's @c NonSpecificEnzymeSearchEngine. Replaces the pre-SNES
     * @c searchDifferentPrecursorRanges flow with a two-phase design:
     *
     *  1. **Byte-count pass**: for every experimental peak, walk the fragment-mz
     *     buckets within fragment tolerance and increment a per-thread byte score
     *     table indexed by @em global mother peptide id. No precursor-mass filter
     *     at this stage — every mother that has a fragment matching any peak is
     *     counted. This is the critical departure from the v1 design, which
     *     pre-filtered mothers by @c mother_mass >= P - tol and admitted the top
     *     half of the index as candidates.
     *
     *  2. **Candidate collection**: for each (precursor charge, isotope error),
     *     walk fragment buckets at the target m/z that a realized sub-peptide's
     *     terminal ion would occupy:
     *       - Single-N mother / b-ion index: target m/z = (M+H)+ − water
     *         (relation @c M_sub+H+ = b_k + water, so the mother's b_k ion falls
     *         at @c M_obs+H+ − water when the realized length matches).
     *       - Single-C mother / y-ion index: target m/z = (M+H)+ directly
     *         (@c M_sub+H+ = y_k exactly).
     *     Every mother with a fragment in the target bin that has the correct
     *     kind (Single-N vs Single-C) and whose byte-score meets
     *     @c fragment:min_matched_ions is emitted as a candidate.
     *
     * This design produces a candidate set sized like a single fragment-bin
     * lookup (dozens, not thousands), matching MetaMorpheus's algorithmic
     * scalability. The byte table is reused across calls via @c thread_local.
     *
     * Only called when @ref isSnesMode is true; otherwise the pre-SNES
     * @c searchDifferentPrecursorRanges path is used.
     *
     * @param[in]  spectrum Experimental spectrum with a single precursor.
     * @param[in]  fasta_entries Source database passed to build(); required
     *             for realizing sub-peptides and applying variable mods
     *             inside the SNES v1.1 subset-enumeration post-pass.
     * @param[out] sms Accumulated candidate matches, ordered by insertion
     *             (caller runs full-score and top-N selection downstream).
     */
    void querySpectrumSNES_(const MSSpectrum& spectrum,
                            const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                            SpectrumMatchesTopN& sms);

    /**
     * @brief queries peaks for a given experimental spectrum with a set range of potential peptides, isotope error and precursor charge. Hits are transferred into a PSM list.
     * Technically an adapter between query(...) and openSearch(...)/searchDifferentPrecursorRanges(...)
     * @param[out] candidates The n best Spectrum matches
     * @param[in] spectrum The queried experimental spectrum
     * @param[in] candidates_range The range of precursors/peptides the peptide could potentially belong to
     * @param[in] isotope_error The applied isotope error
     * @param[in] precursor_charge The applied precursor charge
     */
    void queryPeaks(SpectrumMatchesTopN& candidates,
                   const MSSpectrum& spectrum,
                   const std::pair<size_t, size_t>& candidates_range,
                   const int16_t isotope_error,
                   const uint16_t precursor_charge);
    /**
     * @brief If closed search loops over all isotope errors. For each iteration loop over all peaks with queryPeaks.
     * @brief If open search applies a precursor-mass window
     * @param[in] spectrum experimental query-spectrum
     * @param[in] precursor_mass The mass of the precursor (mz * charge)
     * @param[out] sms The Top m SpectrumMatches
     * @param[in] charge Applied charge
     */
    void searchDifferentPrecursorRanges(const MSSpectrum& spectrum,
                                        float precursor_mass,
                                        SpectrumMatchesTopN& sms,
                                        uint16_t charge);

    /** @brief places the k-largest elements in the front of the input array. Inside of the k-largest elements and outside the elements are not sorted
     *
     */
    void trimHits(SpectrumMatchesTopN& init_hits) const;

    //since we work with TheoreticalSpectrumGenerator, we must transfer some of those member variables
    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;

    // SpectrumGenerator independend member variables
    std::string digestion_enzyme_;
    EnzymaticDigestion::Specificity enzyme_specificity_{EnzymaticDigestion::SPEC_FULL}; ///< 'full' (default), 'semi' (semi-tryptic), or 'none' (e.g. immunopeptidomics)

    size_t missed_cleavages_; ///< number of missed cleavages
    float peptide_min_mass_;
    float peptide_max_mass_;
    size_t peptide_min_length_;
    size_t peptide_max_length_;
  
    StringList modifications_fixed_;    ///< Modification that are one all peptides
    StringList modifications_variable_; ///< Variable Modification -> all possible comibnations are created
    size_t max_variable_mods_per_peptide_;

    // Search Related member variables

    uint16_t min_matched_peaks_;  ///< PSM with less hits are discarded
    int16_t min_isotope_error_;   ///< Minimal possible isotope error
    int16_t max_isotope_error_;   ///< Maximal possible isotope error (both only used for closed search)
    uint16_t min_precursor_charge_; ///< minimal possible precursor charge (usually always 1)
    uint16_t max_precursor_charge_; ///< maximal possible precursor charge
    uint16_t max_fragment_charge_;  ///< The maximal possible charge of the fragments
    uint32_t max_processed_hits_;   ///< The amount of PSM that will be used. the rest is filtered out
    
    /// Instance delegate — same rule, reads the member bounds.
    bool isOpenSearchMode_() const noexcept
    {
      return isOpenSearchMode(precursor_mass_tolerance_lower_,
                              precursor_mass_tolerance_upper_,
                              precursor_mass_tolerance_unit_ppm_);
    }

    /** Compute the signed mass window {lo, hi} around a precursor_mass, converting ppm → Da
     * if the unit is ppm. `lo` is negative (or zero), `hi` is positive (or zero). This is the
     * only place where positive member magnitudes become signed offsets.
     */
    std::pair<float, float> computeMassWindow_(float precursor_mass) const;


  };

}
