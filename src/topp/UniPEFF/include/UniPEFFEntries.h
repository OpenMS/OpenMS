// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include "UniProtXMLFile.h"

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

/**
  @brief The annotations of the PEFF entries that UniPEFF writes for a UniProtKB entry: its canonical entry and the
  entries of its reconstructed isoforms.

  prepareEntry() classifies the features of an entry into annotations of the canonical entry or of the isoform they are
  annotated on; the tool writes the prepared entries.
*/
namespace OpenMS::UniPEFFEntries
{
  // ──────────────────────────────────────────────────────────────────
  // ptmlist.txt — UniProt PTM controlled vocabulary
  // ──────────────────────────────────────────────────────────────────

  struct PtmEntry
  {
    std::string id;               ///< UniProt PTM name (the lookup key from feature/\@description)
    std::string accession;        ///< UniProt internal accession, e.g. "PTM-0253"
    std::string psi_mod;          ///< "MOD:xxxxx" or empty
    std::string unimod;           ///< "N" (numeric Unimod id) or empty (NOT prefixed with "UNIMOD:")
  };

  using PtmMap = std::unordered_map<std::string, PtmEntry>;

  /// Isoforms without a reconstructable sequence ("external": stored as its own
  /// UniProt entry; "not described": sequence unknown), counted per sequence type
  /// so the run reports them in one summary line instead of one line per isoform.
  using SkippedIsoformCounts = std::map<std::string, size_t>;

  // ──────────────────────────────────────────────────────────────────
  // Intermediate annotation structures (built per UniProtEntry)
  // ──────────────────────────────────────────────────────────────────

  /// Sentinel for "no annotation identifier assigned": the tuple is emitted without an id: prefix.
  constexpr uint32_t kNoId = std::numeric_limits<uint32_t>::max();

  struct ModResItem
  {
    int  position{0};        ///< 0 = unknown ('?')
    bool is_halfcys{false};  ///< for stable tie-break on merge
    const PtmEntry* ptm{nullptr};  ///< nullptr = generic (no accession), description in name
    std::string name;        ///< display name (UniProt PTM ID or fallback description)
    uint32_t annotation_id{kNoId};  ///< 0-based; all tuples in Option B, paired half-cystines in default mode
  };

  struct VariantSimpleItem
  {
    int  position{0};
    char new_aa{'?'};
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct VariantComplexItem
  {
    int begin{0};
    int end{0};
    std::string new_seq;
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct ProcessedItem
  {
    int begin{0};
    int end{0};
    std::string cv;          ///< PEFF CV accession, e.g. "PEFF:0001021"
    std::string type_name;   ///< human name, e.g. "signal peptide"
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct DisulfidePairItem
  {
    /// Indices into EntryAnnotations::mods after merge — resolved to annotation_ids during emission.
    size_t idx_a{0};
    size_t idx_b{0};
    bool valid{false};       ///< true when both endpoints come from a &lt;begin&gt;/&lt;end&gt; location; cleared by invalidateOutOfRangeDisulfides()
    uint32_t annotation_id{kNoId};  ///< id of the \\DisulfideBond tuple itself (kNoId = not emitted)
  };

  struct EntryAnnotations
  {
    std::vector<ModResItem> regular_mods;
    std::vector<ModResItem> halfcys;
    std::vector<DisulfidePairItem> disulfides;   ///< indices reference @c halfcys before merge, @c mods after merge
    std::vector<VariantSimpleItem> simple_variants;
    std::vector<VariantComplexItem> complex_variants;
    std::vector<ProcessedItem> processed;

    /// Final merged + sorted list (built by classifyAndMerge).
    std::vector<ModResItem> mods;
  };

  /// One reconstructed isoform, ready for emission as its own PEFF entry.
  struct IsoformEntry
  {
    std::string accession;  ///< isoform accession, e.g. "P02768-2"
    std::string name;       ///< isoform &lt;name&gt;, e.g. "2" or "VEGF121"
    std::string sequence;   ///< canonical sequence with the referenced splice variants applied
    /// features UniProt annotates on this isoform (&lt;location sequence="..."&gt;), plus, with
    /// -map_canonical_features, the canonical ones mapped onto it
    EntryAnnotations annotations;
  };

  struct PreparedEntry
  {
    UniProtEntry source;
    EntryAnnotations annotations;
    size_t isoform_scoped_applied{0};  ///< isoform-scoped features written to their isoform's entry
    size_t isoform_scoped_skipped{0};  ///< isoform-scoped features left out (isoform not written) that would otherwise have become annotations
    /// -map_canonical_features: canonical annotations written to an isoform entry (counted once per isoform)
    size_t canonical_mapped{0};
    /// -map_canonical_features: canonical annotations not written to an isoform entry because they lie in or
    /// cross one of its splice variants, have an unknown position, or contradict an annotation UniProt makes
    /// on that isoform itself (counted once per isoform)
    size_t canonical_unmapped{0};
    std::string prefix;       ///< sp/tr or user override
    std::vector<IsoformEntry> isoforms;  ///< reconstructed alternative products (emitted after the entry)
  };

  /// Counts, over all entries, features that would become annotations of isoform entries: those
  /// written (applied) and those left out (skipped). Used for the isoform-scoped features (skipped:
  /// the isoform is not written: external, not described, not reconstructable, or -omit_isoforms)
  /// and for the canonical features mapped with -map_canonical_features (skipped: in or across a
  /// splice variant, unknown position, or contradicted by the isoform's own annotation).
  struct IsoformScopedTracker
  {
    size_t applied{0};
    size_t applied_entries{0};
    size_t skipped{0};
    size_t skipped_entries{0};

    void add(size_t applied_in_entry, size_t skipped_in_entry)
    {
      applied += applied_in_entry;
      if (applied_in_entry > 0) ++applied_entries;
      skipped += skipped_in_entry;
      if (skipped_in_entry > 0) ++skipped_entries;
    }
  };

  /// PEFF VariantSimple newAA must be an amino-acid letter (or '*').
  bool isResidueCode(char c);

  /// 0 → "?", else decimal.
  std::string positionOrUnknown(int p);

  /// Prepare one entry: reconstruct the emittable isoform sequences ("external" / "not
  /// described" isoforms are counted in @p skipped_isoforms), classify each feature into the
  /// canonical entry or, if its &lt;location sequence="..."&gt; names an emitted isoform, into
  /// that isoform's entry. With @p map_canonical_features, the canonical entry's annotations
  /// are also mapped onto each isoform through its splice variants: features in a deleted or
  /// replaced region, or regions overlapping one, are not mapped (processed regions are clipped
  /// to the residues they keep), and a copy of a feature UniProt annotates on that isoform itself
  /// (for a glycosylation: any on that residue) is not added again. UniProt's own isoform annotation
  /// also wins over a mapped one it contradicts (same-type processed regions that overlap, a disulfide
  /// bond sharing a cysteine, cross-links sharing a residue if one links two residues): the mapped one is
  /// not added and a warning is logged. Then, per entry, merge halfcys, drop disulfide pairs with out-of-range
  /// endpoints and stamp annotation IDs (global in Option B, selective disulfide labels otherwise).
  PreparedEntry prepareEntry(UniProtEntry&& e, const PtmMap& ptms, const std::string& prefix_override,
                             bool option_b, bool record_processing, bool record_aa_mods,
                             bool record_variants, bool record_isoforms, bool map_canonical_features,
                             SkippedIsoformCounts& skipped_isoforms);
} // namespace OpenMS::UniPEFFEntries
