// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
   * @brief The sequence variants of a PEFF database as peptides of a ProSE search (database:peff_variants).
   *
   * A PEFF database read as FASTA (FASTAFile skips the PEFF header lines) keeps the annotations of every entry in its
   * description. For each simple variant (\\VariantSimple: one residue substituted by another) addVariantWindows() appends
   * a @em window to the database: a stretch of the protein around the variant, with the variant residue, long enough that
   * every peptide of the search that the variant changes lies inside it, away from its ends. The window is searched like
   * any other entry, so its peptides are digested with the enzyme, specificity and length limits of the search, and the
   * decoys are generated for it as for the proteins.
   *
   * Each window has a @em region: the residues the variant can change (the variant, the residue before and the two after
   * it, which decide whether the enzyme cleaves next to it). FragmentIndex::setPeptideRegions() with peptideRegions() indexes only the
   * peptides of a window that overlap its region. The others are peptides of the protein itself (already indexed for it)
   * or, at the ends of the window, peptides the protein does not have (the window ends are not protein termini). The
   * decoys of the windows get the corresponding regions (addDecoys()).
   *
   * After the hits are mapped to the entries of the database (PeptideIndexing or the fragment index), mapToParents()
   * maps the hits of the windows back to the proteins they come from: the evidences of a window become evidences of its
   * protein (decoy protein), at the protein positions, and the windows leave the protein hits. A hit whose sequence occurs
   * in a window but in no other entry gets the meta value variantMetaValue() with the variants (e.g. "sp:P04637:R175H";
   * for a decoy, those of the decoy windows, e.g. "DECOY_sp:P04637:R175H"), so that the variant hits can be told apart
   * from the others among the decoys as well (e.g. for a target-decoy competition of their own).
   *
   * How the windows are cut (see addVariantWindows()):
   * - The window reaches at least peptide:max_size + 1 residues beyond the residues the variant can change (in the
   *   target, and with decoys that keep the residues of a cleavage segment in place, in its decoy). A peptide that overlaps
   *   the region therefore never touches an end of the window. With peptide:max_size 0 (no maximum), or 'no cleavage', the
   *   window is the whole protein.
   * - With decoys built per cleavage segment (decoy_method 'reverse' with a cleaving enzyme and a specific digestion, or
   *   'keep_termini') the window starts and ends at cleavage sites of the protein, and one more segment is added at its
   *   end: its decoy then has the residues of the protein's decoy, at the same positions (DecoyGenerator::reversePeptides()
   *   reverses the last segment of a sequence as a whole). Otherwise (the decoy reverses the whole entry) the decoy of the
   *   window is a stretch of the protein's decoy as well, mirrored.
   * - A window does not start with M, unless it starts the protein: a generated decoy keeps an initial M in place when
   *   peptide:clip_nterm_methionine is set, which would shift it against the protein's decoy.
   *
   * Stop codons (\\VariantSimple to a stop codon, `*`, e.g. `Q13*`): the window is the protein up to the residue before
   * the stop, and ends at a protein C-terminus there (its region: the last residue, i.e. the peptides that end at the new
   * C-terminus).
   *
   * Not searched: a variant to the residue it replaces, to a letter without a residue of known mass (X, B, Z, ...), a stop
   * codon at the first residue, a position outside the protein, \\VariantComplex (insertions, deletions) and
   * combinations of variants. getVariantStatus() records every annotated variant with whether it is searched and why
   * not; storeVariantStatus() writes it with the PSMs of each variant.
   */
  class OPENMS_DLLAPI ProSEPEFF
  {
  public:
    /// The search settings that decide where windows are cut
    struct Settings
    {
      std::string enzyme = "Trypsin";                                                     ///< the enzyme of the search
      Size peptide_max_length = 40;                                                       ///< peptide:max_size (0 = no maximum)
      /// Whether the decoys keep the residues of every cleavage segment in place (DecoyGenerator::reversePeptides() with a
      /// specific digestion, or keep_termini); false: the decoy of an entry is the whole entry reversed (non-specific
      /// searches). With 'unspecific cleavage' or 'no cleavage' the decoys are whole entries reversed anyway.
      bool segment_decoys = true;
      /// peptide:clip_nterm_methionine: a generated decoy keeps an initial M in place
      bool clip_nterm_methionine = true;
    };

    /// A window (target or decoy) of the database
    struct Window
    {
      std::string parent;            ///< identifier of the protein (decoy protein) the window comes from
      Size shift = 0;                ///< protein position of window position 0 (window positions + shift = protein positions)
      Size length = 0;               ///< number of residues of the window
      std::pair<Size, Size> region;  ///< [first, second): the window positions a peptide must overlap
      std::string variant;           ///< the variant, e.g. "R175H" (residue, 1-based protein position, variant residue)
      bool decoy = false;            ///< whether the window is the decoy of a window
    };

    /// The kind of an annotated variant
    enum class VariantKind
    {
      SIMPLE,    ///< \\VariantSimple to a residue
      STOP_GAIN, ///< \\VariantSimple to a stop codon (`*`)
      COMPLEX    ///< \\VariantComplex (insertion, deletion, multi-residue substitution)
    };

    /// An annotated variant of the database: whether it is searched, and why not
    struct VariantStatus
    {
      std::string protein;        ///< identifier of the protein
      Size position = 0;          ///< 1-based protein position (a complex variant: its first residue)
      char reference = '\0';      ///< the residue of the protein at the position ('\0': complex, or out of range)
      std::string variant;        ///< e.g. `R175H`, `Q13*`, `120-122|AG` (complex: start-end|replacement)
      VariantKind kind = VariantKind::SIMPLE;
      bool searched = false;      ///< whether a window holds it
      /// why it is not searched (empty if searched): out_of_range, unknown_residue, same_as_reference,
      /// stop_at_first_residue, complex_not_supported
      std::string reason;
      std::string window;         ///< identifier of its window (if searched)
    };

    /// The PSMs of an annotated variant in the hits of one or more runs (see countPSMs())
    struct VariantCounts
    {
      Size psms = 0;                   ///< target PSMs (best hit of a spectrum) with the variant
      Size passing_psms = 0;           ///< ... of them at a q-value <= the threshold (all, without q-values)
      Size passing_runs = 0;           ///< runs with at least one of them
      double best_q = -1.0;            ///< the best q-value of its PSMs (-1: no PSM, or no q-values)
      Size reference_passing_psms = 0; ///< target PSMs without a variant, at q <= the threshold, whose evidence on the
                                       ///< protein covers the variant position with the reference residue there (the
                                       ///< reference sequence is observed)
    };

    /// Name of the meta value of the hits of variant peptides
    static const std::string& variantMetaValue();

    /// Name of a VariantKind as storeVariantStatus() writes it (simple, stop_gain, complex)
    static std::string kindName(VariantKind kind);

    /// A substitution that one of the searched modifications makes identical to the reference (see
    /// referenceEquivalentSubstitutions())
    struct EquivalentSubstitution
    {
      char variant = '\0';      ///< the variant residue, which the modification modifies
      char reference = '\0';    ///< the residue it then equals (same elemental composition)
      std::string modification; ///< the modification, as named in @p modifications
    };

    /**
     * @brief The substitutions reference -> variant that a modification of the variant residue turns back into the
     * reference: the variant residue with the modification has the elemental composition of the reference residue (e.g.
     * E -> Q with Deamidated (Q), D -> N with Deamidated (N)).
     *
     * Such a variant peptide, modified there, is the reference peptide: it matches every spectrum of the reference
     * peptide as well and cannot be told apart from it. Residue modifications only (not terminal ones without residue).
     *
     * @param[in] modifications Modification names as ModifiedPeptideGenerator::getModifications() takes them (e.g. the
     *            fixed and the variable modifications of the search)
     */
    static std::vector<EquivalentSubstitution> referenceEquivalentSubstitutions(const StringList& modifications);

    /// Separator between the identifier of a protein and the variant in the identifier of its windows
    static const std::string& windowSeparator();

    /**
     * @brief Appends a window to @p db for every simple variant of its entries (and remembers the windows).
     *
     * The variants are read from the descriptions of the entries (PEFFFile::parseFASTAEntry()); entries without
     * \\VariantSimple are skipped quickly. A window is named by its protein, windowSeparator() and the variant (e.g.
     * "sp:P04637_VAR_R175H", with a number appended if that identifier exists already). Repeated variants of an entry
     * (same position and residue) give one window.
     *
     * @param[in,out] db The target entries (no decoys): the windows are appended
     * @param[in] settings The settings of the search
     * @return The number of windows added
     */
    Size addVariantWindows(std::vector<FASTAFile::FASTAEntry>& db, const Settings& settings);

    /**
     * @brief Remembers the decoys of the windows: the entries of @p db named @p decoy_prefix followed by the identifier
     * of a window (as ProSE generates them).
     */
    void addDecoys(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& decoy_prefix);

    /// Number of windows (targets and decoys)
    Size size() const { return windows_.size(); }

    /// Whether there are no windows
    bool empty() const { return windows_.empty(); }

    /// The window with identifier @p identifier, or nullptr
    const Window* find(const std::string& identifier) const;

    /// The regions of the windows, for FragmentIndex::setPeptideRegions()
    std::shared_ptr<const FragmentIndex::PeptideRegions> peptideRegions() const { return regions_; }

    /**
     * @brief Maps the hits of windows to their proteins (see class description).
     *
     * For every evidence of a window: dropped if the peptide does not overlap the window's region (a peptide of the
     * protein, which has its own evidence), otherwise moved to the protein (decoy protein) at the protein position (with
     * the residues before and after it in the window). Repeated evidences are merged. 'protein_references' is updated as
     * PeptideIndexing sets it. The protein hits of windows are replaced by hits of their proteins (if not present).
     */
    void mapToParents(std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) const;

    /// Every annotated variant that addVariantWindows() read, in database order, with whether it is searched
    const std::vector<VariantStatus>& getVariantStatus() const { return status_; }

    /**
     * @brief Counts the PSMs of every annotated variant (getVariantStatus(), same order) in the hits of @p runs.
     *
     * The best hit of every spectrum counts (target hits only), after mapToParents(): a hit with variantMetaValue() for
     * the variants it names, a hit without it for the variants whose position one of its evidences covers (on the
     * protein) with the reference residue there (a variant peptide that another entry holds as well carries the variant
     * residue: no reference). The main score is the q-value if the identification's score type is a q-value (rescoring, FDR); otherwise
     * every PSM passes the threshold and best_q stays -1.
     *
     * @param[in] runs The identifications of each run
     * @param[in] q_threshold The q-value threshold of the passing PSMs
     */
    std::vector<VariantCounts> countPSMs(const std::vector<const PeptideIdentificationList*>& runs, double q_threshold) const;

    /**
     * @brief Writes getVariantStatus() with countPSMs() as a tab-separated table.
     *
     * Columns: protein, position, variant, kind, status (searched / not_searched), reason, window, psms,
     * psms_q (passing @p q_threshold), runs_q, best_q (empty without), reference_psms_q.
     *
     * @exception Exception::UnableToCreateFile if @p filename cannot be written
     */
    void storeVariantStatus(const std::string& filename, const std::vector<const PeptideIdentificationList*>& runs,
                            double q_threshold) const;

  private:
    std::unordered_map<std::string, Window> windows_; ///< by identifier
    std::shared_ptr<FragmentIndex::PeptideRegions> regions_ = std::make_shared<FragmentIndex::PeptideRegions>();
    /// The decoy windows to be (addDecoys()), by the identifier of their target window
    std::unordered_map<std::string, Window> pending_decoys_;
    /// getVariantStatus()
    std::vector<VariantStatus> status_;
  };
} // namespace OpenMS
