// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  class PeptideHit;
  class ProteaseDigestion;

  /**
    @brief Knockoff entries for single-residue protein variants, and knockoff q-values for variant PSMs

    A decoy only models random matches. Most false variant calls are not random: a canonical peptide carrying a
    modification of the substitution's mass (N→D and deamidation, S→T and methylation, ...) or picked at the wrong
    isotope fits the variant peptide's spectrum and precursor mass. A <b>knockoff</b> copies a real variant's
    substitution to another position of the same protein that holds the same residue, in a matched sequence context
    and far from every listed variant. It is almost certainly absent from the sample, but a modification or isotope
    error fits it as well as it fits the real variant, so the knockoff PSMs that pass a score estimate the false real
    variant PSMs that pass it.

    <b>Database</b> (buildDatabase()). For a variant at position @em p of a protein with residue @em o, a knockoff
    position @em q
    - holds @em o, and is followed by the same residue as @em p when @em o is N, Q, K or R or either next residue is P
      or G, else by a residue of the same class (hydrophobic AVLIM, aromatic FWY, polar STNQ, basic KRH, acidic DE);
    - is the N of an N-x-S/T sequon exactly when @em p is;
    - lies in another enzymatic peptide (no missed cleavage) of 6 to 50 residues whose length is within 3 of @em p's;
    - is at least 10 residues from every listed variant of the protein.
    Up to @p per_variant positions are drawn, deterministically per variant. A variant entry is named
    <tt>&lt;accession&gt;_VAR_&lt;ref&gt;&lt;pos&gt;&lt;alt&gt;</tt> and a knockoff entry
    <tt>&lt;accession&gt;_KO&lt;i&gt;_&lt;ref&gt;&lt;pos'&gt;&lt;alt&gt;_VAR_&lt;ref&gt;&lt;pos&gt;&lt;alt&gt;</tt> (1-based positions); the
    description carries the variant's prior tier as <tt>tier=&lt;tier&gt;</tt>.

    <b>Knockoff q-values</b> (applyKnockoffQValues()). Rows are the target PSMs whose proteins are all variant or
    knockoff entries. They are split into strata by the substitution's mass (deamidation-like 0.984 Da, methyl-like
    14.016 Da, oxidation-like 15.995 Da, other) and, when variants come in several prior tiers (sample, common,
    curated, catalogue), by tier; a tier with fewer than @p min_tested tested PSMs at the cutoff joins its class's other
    small tiers. Within a stratum, at score @em t:
    FDP(t) = (1/K + Σ 1/k over knockoff rows ≥ t) / #(tested real rows ≥ t),
    with @em k the knockoffs of the variant a knockoff copies and @em K the most knockoffs of any variant. A real row's
    knockoff q-value is the least FDP at or below its score, and its q-value is raised to it. A real row none of whose
    variants has a knockoff gets q = 1: no null covers it. Knockoff rows get q = 1 and are never accepted, and knockoff
    accessions are dropped from every other PSM.

    Peptides that a reference protein holds up to I/L are not variant-specific: index with PeptideIndexer's
    @p IL_equivalent so such PSMs also name the reference protein and are not counted as variant PSMs.

    Ported from bigbio/andes #173, #175 and #180.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI VariantKnockoffs
  {
  public:
    /// A single-residue substitution of a protein
    struct OPENMS_DLLAPI Variant
    {
      std::string accession; ///< protein accession (FASTA identifier, or the middle field of a UniProt 'db|ACC|NAME' identifier)
      Size position = 0;     ///< 0-based position in the protein
      char ref = 'X';        ///< reference residue
      char alt = 'X';        ///< substituted residue
      std::string id;        ///< optional variant id (e.g. a dbSNP or UniProt feature id), written to the description
      std::string tier;      ///< optional prior tier: sample, common, curated or catalogue
    };

    /// What buildDatabase() did
    struct OPENMS_DLLAPI DatabaseStats
    {
      Size variant_entries = 0;           ///< variant entries written
      Size knockoff_entries = 0;          ///< knockoff entries written
      Size variants_without_knockoff = 0; ///< variants for which no position qualified for a knockoff
      Size unknown_protein = 0;           ///< variants of a protein not in the database
      Size reference_mismatch = 0;        ///< variants whose reference residue is not the protein's residue
      Size duplicates = 0;                ///< variants listed more than once
    };

    /// What a variant or knockoff entry is, parsed from its accession and description
    struct OPENMS_DLLAPI EntryInfo
    {
      bool knockoff = false;
      std::string family;   ///< accession of the canonical protein
      std::string variant;  ///< the real variant's substitution, e.g. "D5N" (for a knockoff: the one it copies)
      double delta_mass = 0.0; ///< monoisotopic mass of the substituted minus the reference residue
      std::string tier;     ///< prior tier (may be empty)
    };

    /// The variant and knockoff entries of a database
    struct OPENMS_DLLAPI Index
    {
      std::map<std::string, EntryInfo> entries; ///< by accession
      std::map<std::pair<std::string, std::string>, Size> knockoffs_per_variant; ///< (family, variant) -> knockoffs
      Size max_knockoffs = 0; ///< most knockoffs of any variant (K)
    };

    /// One stratum, as counted at the cutoff
    struct OPENMS_DLLAPI StratumCount
    {
      std::string name;           ///< substitution class, or class/tier, or class/pooled
      Size real_before = 0;       ///< real variant PSMs at q <= cutoff before the knockoff q-values
      Size tested_before = 0;     ///< of them, those of a variant with knockoffs
      double knockoff_weighted = 0.0; ///< knockoff PSMs at q <= cutoff, each counting 1/k
      Size real_after = 0;        ///< real variant PSMs at q <= cutoff after
      /// estimated share of false PSMs among the tested real PSMs at the cutoff (NaN without tested PSMs)
      double falseShare() const;
    };

    /// meta value of a hit: its knockoff q-value (real variant PSMs of a variant with knockoffs)
    static const std::string KNOCKOFF_Q_VALUE;
    /// meta value of a hit: 1 for a knockoff PSM (q-value set to 1)
    static const std::string IS_KNOCKOFF;
    /// meta value of a hit: the stratum its knockoff q-value was estimated in
    static const std::string KNOCKOFF_STRATUM;

    /// Substitution class of a mass difference: "deamidation-like", "methyl-like", "oxidation-like" or "other"
    static std::string substitutionClass(double delta_mass);

    /// Accession of a variant entry (@p position 0-based)
    static std::string variantAccession(const std::string& accession, char ref, Size position, char alt);

    /// Accession of the @p i-th (1-based) knockoff entry, at @p knockoff_position, of a variant at @p position
    static std::string knockoffAccession(const std::string& accession, char ref, Size knockoff_position, Size i, Size position, char alt);

    /**
      @brief Parse a variant or knockoff entry from its accession (and its description, for the tier)
      @return false if @p accession names neither
    */
    static bool parseEntry(const std::string& accession, const std::string& description, EntryInfo& info);

    /**
      @brief Knockoff positions (0-based) for a substitution of residue @p position of @p sequence

      See the class description for the rules. @p listed are the positions of every listed variant of the protein.
      @p digestion splits the protein into enzymatic peptides (its missed cleavages are ignored). Up to
      @p per_variant positions are drawn without replacement by a generator seeded with @p seed.
    */
    static std::vector<Size> knockoffPositions(const std::string& sequence, Size position, const std::vector<Size>& listed,
                                               const ProteaseDigestion& digestion, Size per_variant, UInt64 seed);

    /**
      @brief Variant and knockoff entries for the single-residue variants of the proteins in @p db

      @param db Canonical proteins
      @param variants Variants (their accession must match an entry of @p db)
      @param enzyme Enzyme of the search (knockoffs sit in another enzymatic peptide of a similar length)
      @param per_variant Knockoffs per variant, at most (0: none)
      @param stats Receives counts
      @return The variant entries, followed by the knockoff entries (not the canonical ones)
    */
    static std::vector<FASTAFile::FASTAEntry> buildDatabase(const std::vector<FASTAFile::FASTAEntry>& db,
                                                            const std::vector<Variant>& variants,
                                                            const std::string& enzyme, Size per_variant,
                                                            DatabaseStats& stats);

    /// Index the variant and knockoff entries of a database (decoy entries are not variant entries of a target family)
    static Index indexDatabase(const std::vector<FASTAFile::FASTAEntry>& db);

    /**
      @brief Raise the q-values of variant PSMs to their knockoff q-values

      Must run after the PSM q-values are computed: the main score of each identification must be a q-value. Decoy
      hits are ignored. Does nothing if @p index has no knockoff entry. See the class description.

      @param ids Identifications (main score: q-value); re-sorted afterwards
      @param index Variant and knockoff entries of the searched database
      @param rank_of Score of a hit by which rows are ranked (higher is better)
      @param cutoff q-value at which the strata are counted and small tiers are pooled
      @param min_tested Least tested PSMs at the cutoff for a tier to get its own estimate
      @param use_all_hits Use all hits, not only the best hit of each identification
      @return Per stratum, in a fixed order
    */
    static std::vector<StratumCount> applyKnockoffQValues(PeptideIdentificationList& ids, const Index& index,
                                                          const std::function<double(const PeptideHit&)>& rank_of,
                                                          double cutoff, Size min_tested = 50, bool use_all_hits = false);
  };
} // namespace OpenMS
