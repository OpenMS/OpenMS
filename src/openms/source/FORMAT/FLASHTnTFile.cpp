// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <algorithm>
#include <exception>
#include <sstream>
#include <unordered_map>

namespace OpenMS
{
/**
  @brief FLASHTnT output *.tsv file format
   @ingroup FileIO
**/
void FLASHTnTFile::writeTagHeader(std::fstream& fs)
{
  fs << "TagIndex\tScan\tRetentionTime\tProteoformIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPosition\tDeltaMass\tL"
        "ength\tDeNovoScore\tMasses\tMassScores\n";
}

/// write header line for PrSM file
void FLASHTnTFile::writePrSMHeader(std::fstream& fs, bool chimeric, bool targets)
{
  fs << "PrSMIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tProteoformLevel\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue"
        "\tProteoformIndex\tCandidateSource\tFDRClass\tModVariable\tDecoyPrecursorWeight\tProteinLevelQvalue"
        "\tExtensionScore\tRawScore\tRawMatchedFragments\tRawFragmentHypotheses\tRawCoverage(%)"
     << (chimeric ? "\tChimericRound\tMatchedMasses\tRemovedMasses\tResidualNullScore\tResidualScoreRatio\tUniqueMasses\tUniqueScore" : "")
     << (targets ? "\tTargetProtein\tDecoyMargin" : "") << "\n";
}

/// write header line for Proteoform file
void FLASHTnTFile::writeProHeader(std::fstream& fs, bool chimeric, bool targets)
{
  fs << "ProteoformIndex\tPrSMIndices\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tProteoformLevel\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tProteoformLevelQvalue\tCandidateSource\tFDRClass\tModVariable\tProteinLevelQvalue"
        "\tExtensionScore\tRawScore\tRawMatchedFragments\tRawFragmentHypotheses\tRawCoverage(%)"
     << (chimeric ? "\tChimericRound\tMatchedMasses\tRemovedMasses\tResidualNullScore\tResidualScoreRatio\tUniqueMasses\tUniqueScore" : "") << "\tCombinedPrSMs\tCombinedCoverage(%)\tCombinedCoverageChance(%)" << (targets ? "\tTargetProtein\tDecoyMargin" : "") << "\n";
}

/// write the features in regular file output
void FLASHTnTFile::writeProteinHeader(std::fstream& fs)
{
  fs << "Rank\tScan\tProteinAccession\tProteinDescription\tProteinLength\tScore\tMatchedTagCount\tMatchedPositionCount\tTagIndices\n";
}

void FLASHTnTFile::writeProteins(const std::vector<ProteinHit>& hits, std::fstream& fs)
{
  // one line per candidate protein of every spectrum (hundreds of thousands): formatted in parallel blocks, written in order; the meta
  // values by index, because a lookup by name takes the registry's critical section
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_tag_indices = registry.registerName("TagIndices"), mi_rank = registry.registerName("ProteinRank"), mi_scan = registry.registerName("Scan"),
             mi_matched_aa = registry.registerName("MatchedAA"), mi_description = registry.registerName("Description"),
             mi_protein_length = registry.registerName("ProteinLength");
  const SignedSize block_size = 2048;
  const SignedSize block_count = ((SignedSize)hits.size() + block_size - 1) / block_size;
  std::vector<std::string> lines(block_count);
  std::vector<std::exception_ptr> errors(block_count);
#pragma omp parallel for schedule(dynamic)
  for (SignedSize b = 0; b < block_count; b++)
  {
    try
    {
      std::stringstream ss; // the default formatting, as that of the file stream
      for (SignedSize i = b * block_size; i < std::min((SignedSize)hits.size(), (b + 1) * block_size); i++)
      {
        const auto& hit = hits[i];
        const std::vector<int> tag_indices = hit.getMetaValue(mi_tag_indices).toIntList();
        ss << (int)hit.getMetaValue(mi_rank) << "\t" << (int)hit.getMetaValue(mi_scan) << "\t" << hit.getAccession() << "\t" << hit.getMetaValue(mi_description).toString()
           << "\t" << (hit.metaValueExists(mi_protein_length) ? (Size)(int)hit.getMetaValue(mi_protein_length) : hit.getSequence().size()) << "\t"
           << hit.getScore() << "\t" << tag_indices.size() << "\t" << (int)hit.getMetaValue(mi_matched_aa) << "\t";
        for (int ti : tag_indices) { ss << ti << ";"; }
        ss << "\n";
      }
      lines[b] = ss.str();
    }
    catch (...)
    {
      errors[b] = std::current_exception();
    }
  }
  // an error is rethrown before anything is written, as by the serial loop
  for (const auto& e : errors)
  {
    if (e) std::rethrow_exception(e);
  }
  for (const auto& l : lines) fs << l;
}

void FLASHTnTFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
{
  const auto& tags = tnt.getTags();
  const auto& proteoforms = tnt.getProteoforms();
  const std::vector<int> no_hits;
  // meta values by index: a lookup by name takes the registry's critical section, in every thread below
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_description = registry.registerName("Description"), mi_index = registry.registerName("Index"),
             mi_start = registry.registerName("StartPosition"), mi_end = registry.registerName("EndPosition"),
             mi_mod_variable = registry.registerName("ModificationVariable"), mi_modifications = registry.registerName("Modifications"),
             mi_mod_starts = registry.registerName("ModificationStarts"), mi_entry_start = registry.registerName("EntryStart"),
             mi_entry_end = registry.registerName("EntryEnd");

  // the lines of the tags with PrSMs first, then those of the others, each in tag order; formatted in parallel by blocks of tags
  const int tag_count = (int)tags.size();
  const int block_size = 4096;
  const int block_count = (tag_count + block_size - 1) / block_size;
  std::vector<std::string> with_hits(block_count), without_hits(block_count);
  std::vector<std::exception_ptr> errors(block_count);
#pragma omp parallel for schedule(dynamic)
  for (int b = 0; b < block_count; b++)
  {
    try
    {
      std::stringstream ss_with, ss_without;
      for (int t = b * block_size; t < std::min(tag_count, (b + 1) * block_size); t++)
      {
        const auto& tag = tags[t];
        const std::vector<int>* matched = tag.getIndex() < 0 ? nullptr : tnt.getProteoformIndicesMatchedBy(tag.getIndex());
        const auto& hit_indices = matched == nullptr ? no_hits : *matched;
        std::stringstream& ss = hit_indices.empty() ? ss_without : ss_with;

        std::string acc = "";
        std::string description = "";
        std::string hitindices = "";
        std::string positions = "";
        std::string delta_masses = "";
        for (const int hit_index : hit_indices)
        {
          const auto& hit = proteoforms[hit_index];
          if (! hitindices.empty()) {
            acc += ";";
            description += ";";
            hitindices += ";";
            positions += ";";
            delta_masses += ";";
          }

          acc += hit.getAccession();

          std::string proteindescription = hit.getMetaValue(mi_description).toString(); // as getDescription()
          if (proteindescription.empty()) { proteindescription = " "; }
          description += proteindescription;
          hitindices += StringUtils::toStr(hit.getMetaValue(mi_index));

          auto pos = std::vector<int>();
          auto masses = std::vector<double>();
          auto pos_in_truncated = std::vector<int>();
          auto masses_in_truncated = std::vector<double>();

          int protein_start_position = hit.getMetaValue(mi_start);
          int protein_end_position = hit.getMetaValue(mi_end);
          protein_start_position--;
          std::string seq = hit.getSequence();
          // the proteoform's residues [first, last): an open terminus bounded by the database entry's processed form (EntryStart, EntryEnd
          // of a PEFF entry; else the protein's ends)
          const int first = protein_start_position >= 0 ? protein_start_position
                                                        : (hit.metaValueExists(mi_entry_start) ? std::max(0, (int)hit.getMetaValue(mi_entry_start) - 1) : 0);
          const int last = protein_end_position >= 0 ? protein_end_position
                                                     : (hit.metaValueExists(mi_entry_end) ? (int)hit.getMetaValue(mi_entry_end) : (int)seq.size());
          // the fixed modifications of a PEFF entry (ModificationVariable 2) are part of the flanking residue masses: on the protein, and
          // on the proteoform's residues
          std::vector<std::pair<int, double>> fixed, fixed_in_truncated;
          if (hit.metaValueExists(mi_mod_variable))
          {
            const std::vector<int> variable = hit.getMetaValue(mi_mod_variable).toIntList();
            if (std::find(variable.begin(), variable.end(), 2) != variable.end())
            {
              const std::vector<double> fixed_masses = hit.getMetaValue(mi_modifications);
              const std::vector<int> fixed_starts = hit.getMetaValue(mi_mod_starts).toIntList();
              for (Size m = 0; m < variable.size() && m < fixed_masses.size() && m < fixed_starts.size(); m++)
              {
                if (variable[m] != 2) continue;
                fixed.emplace_back(fixed_starts[m], fixed_masses[m]);
                if (fixed_starts[m] >= first && fixed_starts[m] < last) fixed_in_truncated.emplace_back(fixed_starts[m] - first, fixed_masses[m]);
              }
            }
          }
          FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos, masses, -1, seq, tag, fixed.empty() ? nullptr : &fixed);

          if (last < (int)seq.size()) seq = seq.substr(0, std::max(0, last));
          if (first > 0) seq = seq.substr(std::min((int)seq.size(), first));
          FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos_in_truncated, masses_in_truncated, flanking_mass_tol, seq, tag,
                                                                         fixed.empty() ? nullptr : &fixed_in_truncated);
          if (!pos_in_truncated.empty())
          {
            for (int i = 0; i < pos.size(); i++)
            {
              if (pos[i] - first != pos_in_truncated[0]) continue;
              positions += std::to_string(pos[i] + 1);
              delta_masses += std::to_string(masses[i]);
              break;
            }
          }
        }

        ss << tag.getIndex() << "\t" << tag.getScan() << "\t" << tag.getRetentionTime() << "\t" << hitindices << "\t" << acc << "\t" << description
           << "\t" << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
           << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

        for (const auto& mz : tag.getMzs())
        {
          ss << std::to_string(mz) << ",";
        }
        ss << "\t";
        for (size_t i = 0; i < tag.getMzs().size(); i++) // one score per mass (a gapped tag spans more residues than masses)
        {
          ss << std::to_string(tag.getScore(i)) << ",";
        }
        ss << "\n";
      }
      with_hits[b] = ss_with.str();
      without_hits[b] = ss_without.str();
    }
    catch (...)
    {
      errors[b] = std::current_exception();
    }
  }
  // an error is rethrown before anything is written, as by the serial loop
  for (const auto& e : errors)
  {
    if (e) std::rethrow_exception(e);
  }
  for (const auto& lines : with_hits) fs << lines;
  for (const auto& lines : without_hits) fs << lines;
}

namespace
{
  /// the registry index of a meta value name of the PrSM output, registered once: a lookup by name takes the registry's critical section,
  /// which serialised the writers (the map is read-only after its construction, so threads share it without a lock)
  OpenMS::UInt prsmKey(const char* name)
  {
    static const std::unordered_map<std::string, OpenMS::UInt> keys = [] {
      auto& registry = OpenMS::MetaInfoInterface::metaRegistry();
      std::unordered_map<std::string, OpenMS::UInt> k;
      for (const char* n : {"TargetProtein", "DecoyMargin", "ChimericRound", "CombinedCoverage", "CombinedCoverageChance", "CombinedPrSMs", "DecoyPrecursor", "Description", "EndPosition", "ExtensionScore", "FDRClass", "GappedTag", "GivenMass", "Index", "Mass", "MatchedAA", "MatchedMassIndices", "ModificationACCs", "ModificationEnds", "ModificationIDs", "ModificationStarts", "ModificationVariable", "Modifications", "NumMass", "PrecursorLadder", "PrecursorSNR", "PrecursorScore", "ProteinQvalue", "ProteoformIndex", "ProteoformLevel", "ProteoformMassByFragmentMass", "RT", "RawCoverage", "RawFragmentHypotheses", "RawMatchedFragments", "RawScore", "RemovedMasses", "ResidualNullScore", "ResidualScoreRatio", "Scan", "SlidingLadder", "StartPosition", "TagIndices", "TerminalLadder", "TheoreticalMass", "UniqueMasses", "UniqueScore", "proqvalue", "qvalue"})
        k.emplace(n, registry.registerName(n));
      return k;
    }();
    return keys.at(name);
  }

  /// the proteoform mass of the output: of a hit with a mass its own (TheoreticalMass: residues, fixed and variable modifications at their
  /// masses, blind shifts as measured), else Mass (-1: none). Mass itself is the measured one the search compares.
  OpenMS::DataValue reportedMass(const OpenMS::ProteinHit& hit)
  {
    const OpenMS::DataValue& mass = hit.getMetaValue(prsmKey("Mass"));
    if ((double)mass <= 0 || ! hit.metaValueExists(prsmKey("TheoreticalMass"))) return mass;
    const OpenMS::DataValue& theoretical = hit.getMetaValue(prsmKey("TheoreticalMass"));
    return (double)theoretical > 0 ? theoretical : mass;
  }
} // namespace

namespace
{
  /// where the candidate protein of a hit came from
  std::string candidateSource(const OpenMS::ProteinHit& hit)
  {
    if (hit.metaValueExists(prsmKey("PrecursorLadder"))) return "precursor ladder";
    if (hit.metaValueExists(prsmKey("TerminalLadder"))) return "terminal ladder";
    if (hit.metaValueExists(prsmKey("SlidingLadder"))) return "sliding ladder";
    if (hit.metaValueExists(prsmKey("GappedTag"))) return "gapped tag";
    return "sequence tag";
  }

  /// 1 for a variable (searched) modification, 0 for a blind shift, 2 for a fixed one of the database entry (PEFF), per modification
  std::string modVariable(const OpenMS::ProteinHit& hit)
  {
    std::string out;
    if (! hit.metaValueExists(prsmKey("ModificationVariable"))) return out;
    for (int v : hit.getMetaValue(prsmKey("ModificationVariable")).toIntList())
    {
      if (! out.empty()) out += ";";
      out += v == 2 ? "2" : (v != 0 ? "1" : "0");
    }
    return out;
  }

  /// the residues [first, end) of the hit's sequence that the ProteinSequence column shows: the proteoform, an open terminus bounded by
  /// the database entry's processed form (EntryStart, EntryEnd of a PEFF entry; else the protein's ends)
  std::pair<int, int> proteinSequenceRange(const OpenMS::ProteinHit& hit)
  {
    const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
    const int size = (int)hit.getSequence().size();
    const int entry_first = hit.metaValueExists("EntryStart") ? (int)hit.getMetaValue("EntryStart") - 1 : 0;
    const int entry_end = hit.metaValueExists("EntryEnd") ? (int)hit.getMetaValue("EntryEnd") : size;
    const int first = start < 0 ? entry_first : (start - 1);
    const int last = end < 0 ? entry_end : end;
    return {first, last};
  }

  /// the raw-evidence columns: extension score, raw score, matched fragments, fragment hypotheses, coverage (-1 without a raw file)
  std::string rawColumns(const OpenMS::ProteinHit& hit)
  {
    std::string out = "\t";
    out += StringUtils::toStr(hit.metaValueExists(prsmKey("ExtensionScore")) ? (double)hit.getMetaValue(prsmKey("ExtensionScore")) : hit.getScore());
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("RawScore")) ? (double)hit.getMetaValue(prsmKey("RawScore")) : -1.0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("RawMatchedFragments")) ? (int)hit.getMetaValue(prsmKey("RawMatchedFragments")) : -1);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("RawFragmentHypotheses")) ? (int)hit.getMetaValue(prsmKey("RawFragmentHypotheses")) : -1);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("RawCoverage")) ? (double)hit.getMetaValue(prsmKey("RawCoverage")) : -1.0);
    return out;
  }

  /// the combined-evidence columns of a proteoform (FLASHTnTAlgorithm::setCombinedCoverage_): its PrSMs combined, their coverage and its
  /// chance (-1 without a raw file)
  std::string combinedColumns(const OpenMS::ProteinHit& hit)
  {
    std::string out = "\t";
    out += StringUtils::toStr(hit.metaValueExists(prsmKey("CombinedPrSMs")) ? (int)hit.getMetaValue(prsmKey("CombinedPrSMs")) : -1);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("CombinedCoverage")) ? (double)hit.getMetaValue(prsmKey("CombinedCoverage")) : -1.0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("CombinedCoverageChance")) ? (double)hit.getMetaValue(prsmKey("CombinedCoverageChance")) : -1.0);
    return out;
  }

  /// the chimeric columns (FLASHTnT -chimeric_rounds): round, matched masses, removed masses, residual null score, score ratio, and of
  /// native mode the masses of a PrSM no better PrSM of its spectrum matched and their chance score (-1: the first PrSM, or not native)
  std::string chimericColumns(const OpenMS::ProteinHit& hit)
  {
    std::string out = "\t";
    out += StringUtils::toStr(hit.metaValueExists(prsmKey("ChimericRound")) ? (int)hit.getMetaValue(prsmKey("ChimericRound")) : 0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("MatchedMassIndices")) ? (int)hit.getMetaValue(prsmKey("MatchedMassIndices")).toIntList().size() : -1);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("RemovedMasses")) ? (int)hit.getMetaValue(prsmKey("RemovedMasses")) : 0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("ResidualNullScore")) ? (double)hit.getMetaValue(prsmKey("ResidualNullScore")) : -1.0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("ResidualScoreRatio")) ? (double)hit.getMetaValue(prsmKey("ResidualScoreRatio")) : -1.0);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("UniqueMasses")) ? (int)hit.getMetaValue(prsmKey("UniqueMasses")) : -1);
    out += "\t" + StringUtils::toStr(hit.metaValueExists(prsmKey("UniqueScore")) ? (double)hit.getMetaValue(prsmKey("UniqueScore")) : -1.0);
    return out;
  }

  /// the target columns (FLASHTnT -target): TargetProtein (0: not given with -target; 1: given, and reported by the search; 2: given, and
  /// characterized on request only) and DecoyMargin (the candidate score minus that of the best decoy candidate of its source in the
  /// spectrum; empty: not measured)
  std::string targetColumns(const OpenMS::ProteinHit& hit)
  {
    std::string out = "\t";
    out += StringUtils::toStr(hit.metaValueExists(prsmKey("TargetProtein")) ? (int)hit.getMetaValue(prsmKey("TargetProtein")) : 0);
    out += "\t";
    if (hit.metaValueExists(prsmKey("DecoyMargin"))) out += StringUtils::number((double)hit.getMetaValue(prsmKey("DecoyMargin")), 2);
    return out;
  }
} // namespace

namespace
{
  std::string global_fixed_proforma; // FLASHTnTFile::setGlobalFixedModifications
}

void FLASHTnTFile::setGlobalFixedModifications(const std::string& proforma_prefix)
{
  global_fixed_proforma = proforma_prefix;
}

std::string FLASHTnTFile::getProFormaString(const ProteinHit& hit)
{
  const int start = hit.getMetaValue(prsmKey("StartPosition")), end = hit.getMetaValue(prsmKey("EndPosition"));
  const std::string& sequence = hit.getSequence();
  if (start <= 0 || end <= 0 || end < start || end > (int)sequence.size()) return "";
  const int offset = start - 1, length = end - offset;
  const std::string proteoform = sequence.substr(offset, length);
  std::vector<FLASHHelperClasses::ProFormaModification> pf_mods; // the modifications written (ProForma::toString)
  auto add = [&pf_mods](int s, int e, double mass, const std::string& tag) {
    FLASHHelperClasses::ProFormaModification mod;
    mod.start = s;
    mod.end = e;
    mod.mass = mass;
    mod.tag = tag;
    pf_mods.push_back(mod);
  };

  struct Entry
  {
    int start, end;
    double mass;
    std::string names;
    bool variable;
    std::string fixed_accession; // a fixed modification of the database entry (ModificationVariable 2): its accession; else empty
  };
  std::vector<Entry> entries;
  if (hit.metaValueExists(prsmKey("Modifications")))
  {
    const std::vector<double> masses = hit.getMetaValue(prsmKey("Modifications"));
    const std::vector<int> starts = hit.getMetaValue(prsmKey("ModificationStarts"));
    const std::vector<int> ends = hit.getMetaValue(prsmKey("ModificationEnds"));
    const std::vector<std::string> ids = hit.metaValueExists(prsmKey("ModificationIDs")) ? hit.getMetaValue(prsmKey("ModificationIDs")).toStringList() : std::vector<std::string>();
    const std::vector<std::string> accs = hit.metaValueExists(prsmKey("ModificationACCs")) ? hit.getMetaValue(prsmKey("ModificationACCs")).toStringList() : std::vector<std::string>();
    const std::vector<int> variable = hit.metaValueExists(prsmKey("ModificationVariable")) ? (std::vector<int>)hit.getMetaValue(prsmKey("ModificationVariable")).toIntList() : std::vector<int>();
    for (Size i = 0; i < masses.size() && i < starts.size() && i < ends.size(); i++)
    {
      const int s = std::max(0, std::min(length - 1, std::min(starts[i], ends[i]) - offset));
      const int e = std::max(0, std::min(length - 1, std::max(starts[i], ends[i]) - offset));
      const bool fixed = i < variable.size() && variable[i] == 2;
      entries.push_back({s, e, masses[i], i < ids.size() ? ids[i] : std::string(), i < variable.size() && variable[i] == 1,
                         fixed ? (i < accs.size() && ! accs[i].empty() ? accs[i] : std::string("-")) : std::string()});
    }
  }
  std::stable_sort(entries.begin(), entries.end(), [](const Entry& l, const Entry& r) { return l.start < r.start; });
  // ProForma ranges are the same or disjoint: two ranges that overlap otherwise become one, with their summed mass
  std::vector<Entry> merged;
  for (const auto& entry : entries)
  {
    bool combined = false;
    for (auto& m : merged)
    {
      if (entry.end == entry.start || m.end == m.start || (m.start == entry.start && m.end == entry.end)) continue;
      if (m.start <= entry.end && entry.start <= m.end)
      {
        m.start = std::min(m.start, entry.start);
        m.end = std::max(m.end, entry.end);
        m.mass += entry.mass;
        m.names = "";
        m.variable = false;
        m.fixed_accession = "";
        combined = true;
        break;
      }
    }
    if (! combined) merged.push_back(entry);
  }
  for (const auto& m : merged)
  {
    // a fixed modification of the database entry: its mass with its accession as an alternative (e.g. K[+42.010565|MOD:00064];
    // ProteoformDatabase::proformaTag)
    if (! m.fixed_accession.empty())
    {
      add(m.start, m.end, m.mass, ProteoformDatabase::proformaTag(m.fixed_accession, m.mass));
      continue;
    }
    // the names of a variable modification ("A," or "A+B,"; alternatives "A,B," leave the identity open)
    std::vector<std::pair<std::string, double>> named;
    std::string names = m.names;
    while (StringUtils::hasSuffix(names, ",")) names.pop_back();
    if (m.variable && ! names.empty() && ! StringUtils::has(names, ','))
    {
      std::vector<std::string> parts;
      StringUtils::split(names, '+', parts);
      double sum = 0;
      for (const auto& part : parts)
      {
        try
        {
          const double mass = FLASHHelperClasses::getProFormaTagMass(part, proteoform[m.start]);
          named.emplace_back(part, mass);
          sum += mass;
        }
        catch (const Exception::ParseError&)
        {
          named.clear();
          break;
        }
      }
      if (std::abs(sum - m.mass) > 0.05) named.clear();
    }
    if (named.empty()) add(m.start, m.end, m.mass, "");
    for (const auto& [name, mass] : named) add(m.start, m.end, mass, name);
  }
  try
  {
    return global_fixed_proforma + ProForma::toString(FLASHHelperClasses::toPeptidoform(proteoform, pf_mods));
  }
  catch (const Exception::ParseError&)
  {
    return ""; // a residue ProForma has no letter for
  }
  catch (const Exception::InvalidValue&)
  {
    return ""; // ranges that overlap after merging
  }
}


void OpenMS::FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs, bool chimeric, bool targets)
{
  // the rows formatted in parallel, each into its own string, written in input order; an error is rethrown before anything is written, as
  // by the serial loop
  prsmKey("Index"); // the keys registered before the threads
  std::vector<std::string> rows(hits.size());
  std::vector<std::exception_ptr> errors(hits.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (SignedSize h = 0; h < (SignedSize)hits.size(); h++)
  {
    const auto& hit = hits[h];
    if (! hit.metaValueExists(prsmKey("Index"))) continue;
    try
    {
    std::stringstream ss;
    std::string tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists(prsmKey("TagIndices"))) indices = (std::vector<int>)hit.getMetaValue(prsmKey("TagIndices")).toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue(prsmKey("Modifications"));
    std::vector<int> mod_starts = hit.getMetaValue(prsmKey("ModificationStarts"));
    std::vector<int> mod_ends = hit.getMetaValue(prsmKey("ModificationEnds"));
    std::vector<std::string> mod_ids = hit.getMetaValue(prsmKey("ModificationIDs"));
    std::vector<std::string> mod_accs = hit.getMetaValue(prsmKey("ModificationACCs"));

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }

    int start = hit.getMetaValue(prsmKey("StartPosition"));
    int end = hit.getMetaValue(prsmKey("EndPosition"));

    const auto [start_in_seq, end_in_seq] = proteinSequenceRange(hit);

    // Use ProForma for sequence generation
    std::string proformaStr = getProFormaString(hit);
    ss << hit.getMetaValue(prsmKey("Index")) << "\t" << hit.getMetaValue(prsmKey("Scan")) << "\t" << hit.getMetaValue(prsmKey("RT")) << "\t" << hit.getMetaValue(prsmKey("NumMass"))
       << "\t" << hit.getAccession() << "\t" << hit.getMetaValue(prsmKey("Description")).toString() << "\t" << hit.getMetaValue(prsmKey("GivenMass"))  << "\t" << reportedMass(hit) << "\t" << ((int)hit.getMetaValue(prsmKey("ProteoformMassByFragmentMass")) > 0) << "\t" <<  hit.getSequence() << "\t"
       << StringUtils::substr(hit.getSequence(), start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << (hit.metaValueExists(prsmKey("ProteoformLevel")) ? hit.getMetaValue(prsmKey("ProteoformLevel")).toString() : std::string()) << "\t" << hit.getMetaValue(prsmKey("MatchedAA")) << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue(prsmKey("PrecursorScore"))  << "\t" << hit.getMetaValue(prsmKey("PrecursorSNR"))  << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists(prsmKey("qvalue")) ? (double)hit.getMetaValue(prsmKey("qvalue")) : -1)) << "\t"
       << std::to_string((hit.metaValueExists(prsmKey("proqvalue")) ? (double)hit.getMetaValue(prsmKey("proqvalue")) : -1)) << "\t"
       << (hit.metaValueExists(prsmKey("ProteoformIndex")) ? (int)hit.getMetaValue(prsmKey("ProteoformIndex")) : -1) << "\t" << candidateSource(hit) << "\t"
       << (hit.metaValueExists(prsmKey("FDRClass")) ? hit.getMetaValue(prsmKey("FDRClass")).toString() : std::string()) << "\t" << modVariable(hit) << "\t"
       << (hit.metaValueExists(prsmKey("DecoyPrecursor")) ? (double)hit.getMetaValue(prsmKey("DecoyPrecursor")) : 0.0) << "\t"
       << std::to_string(hit.metaValueExists(prsmKey("ProteinQvalue")) ? (double)hit.getMetaValue(prsmKey("ProteinQvalue")) : -1) << rawColumns(hit)
       << (chimeric ? chimericColumns(hit) : std::string()) << (targets ? targetColumns(hit) : std::string()) << "\n";
    rows[h] = ss.str();
    }
    catch (...)
    {
      errors[h] = std::current_exception();
    }
  }
  for (const auto& e : errors)
  {
    if (e) std::rethrow_exception(e);
  }
  for (const auto& row : rows) fs << row;
}

void OpenMS::FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr, bool chimeric, bool targets)
{
  std::stringstream ss;
  int pro_index = 0;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    if ((int)hit.getMetaValue("Representative") == 0) continue;
    if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
    std::string prsmindices = "", tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    //std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> tag_indices;
    if (hit.metaValueExists("TagIndices")) tag_indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : tag_indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<std::string> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<std::string> mod_accs = hit.getMetaValue("ModificationACCs");

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }
    if (hit.metaValueExists("PrSMIndices"))
    {
      auto prsm_indices = hit.getMetaValue("PrSMIndices").toIntList();
      for (const int& index : prsm_indices)
      {
        if (! prsmindices.empty()) prsmindices += ";";
        prsmindices += std::to_string(index);
      }
    }
    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    const auto [start_in_seq, end_in_seq] = proteinSequenceRange(hit);

    // Use ProForma to generate ProForma string
    std::string proformaStr = getProFormaString(hit);

    const int this_index = pro_index++;
    ss << (hit.metaValueExists("ProteoformIndex") ? (int)hit.getMetaValue("ProteoformIndex") : this_index) << "\t" << prsmindices << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass") << "\t" << reportedMass(hit) << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" << hit.getSequence() << "\t"
       << StringUtils::substr(hit.getSequence(), start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << (hit.metaValueExists("ProteoformLevel") ? hit.getMetaValue("ProteoformLevel").toString() : std::string()) << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore") << "\t" << hit.getMetaValue("PrecursorSNR") << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\t" << candidateSource(hit) << "\t"
       << (hit.metaValueExists("ProteoformFDRClass") ? hit.getMetaValue("ProteoformFDRClass").toString()
          : hit.metaValueExists("FDRClass") ? hit.getMetaValue("FDRClass").toString() : std::string()) << "\t" << modVariable(hit) << "\t"
       << std::to_string(hit.metaValueExists("ProteinQvalue") ? (double)hit.getMetaValue("ProteinQvalue") : -1) << rawColumns(hit)
       << (chimeric ? chimericColumns(hit) : std::string()) << combinedColumns(hit) << (targets ? targetColumns(hit) : std::string()) << "\n";
  }
  fs << ss.str();
}

} // namespace OpenMS
// namespace OpenMS