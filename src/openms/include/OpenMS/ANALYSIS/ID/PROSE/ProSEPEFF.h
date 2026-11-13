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
   * in a window but in no other entry gets the meta value variantMetaValue() with the variants (e.g. "sp:P04637:R175H").
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
   * Simple variants only: a variant to the residue it replaces, to a stop codon or to a letter without a residue of known
   * mass (X, B, Z, ...) is skipped, and so are \\VariantComplex (insertions, deletions) and combinations of variants.
   */
  class OPENMS_DLLAPI ProSEPEFF
  {
  public:
    /// The search settings that decide where windows are cut
    struct Settings
    {
      std::string enzyme = "Trypsin";                                                     ///< the enzyme of the search
      EnzymaticDigestion::Specificity specificity = EnzymaticDigestion::SPEC_FULL;        ///< its specificity
      Size peptide_max_length = 40;                                                       ///< peptide:max_size (0 = no maximum)
      /// Whether the decoys keep the residues of every cleavage segment in place (DecoyGenerator::reversePeptides(), or
      /// keep_termini); false: the decoy of an entry is the whole entry reversed (non-specific searches, 'unspecific
      /// cleavage')
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

    /// Name of the meta value of the hits of variant peptides
    static const std::string& variantMetaValue();

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

  private:
    std::unordered_map<std::string, Window> windows_; ///< by identifier
    std::shared_ptr<FragmentIndex::PeptideRegions> regions_ = std::make_shared<FragmentIndex::PeptideRegions>();
    /// The decoy windows to be (addDecoys()), by the identifier of their target window
    std::unordered_map<std::string, Window> pending_decoys_;
  };
} // namespace OpenMS
