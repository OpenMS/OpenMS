// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/CHEMISTRY/ProForma.h> // added for ProForma
#include <OpenMS/FORMAT/FLASHTnTFile.h>

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
void FLASHTnTFile::writePrSMHeader(std::fstream& fs)
{
  fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue\n";
}

/// write header line for Proteoform file
void FLASHTnTFile::writeProHeader(std::fstream& fs)
{
  fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tProteoformLevelQvalue\n";
}

/// write the features in regular file output
void FLASHTnTFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
{
  std::stringstream ss;
  auto tags = std::vector<FLASHHelperClasses::Tag>();
  tnt.getTags(tags);

  for (int c = 0; c < 2; c++)
  {
    for (const auto& tag : tags)
    {
      auto hits = std::vector<ProteinHit>();
      tnt.getProteoformHitsMatchedBy(tag, hits);
      if (c == 0 && hits.empty()) continue;
      if (c == 1 && ! hits.empty()) continue;

      String acc = "";
      String description = "";
      String hitindices = "";
      String positions = "";
      String delta_masses = "";
      for (const auto& hit : hits)
      {
        if (! acc.empty()) acc += ";";
        if (! description.empty()) description += ";";
        if (! hitindices.empty()) hitindices += ";";
        if (! positions.empty()) positions += ";";
        if (! delta_masses.empty()) delta_masses += ";";

        acc += hit.getAccession();

        String proteindescription = hit.getDescription();
        if (proteindescription.empty()) { proteindescription = " "; }
        description += proteindescription;
        hitindices += (String)hit.getMetaValue("Index");

        auto pos = std::vector<int>();
        auto masses = std::vector<double>();
        auto pos_in_truncated = std::vector<int>();
        auto masses_in_truncated = std::vector<double>();

        int protein_start_position = hit.getMetaValue("StartPosition");
        int protein_end_position = hit.getMetaValue("EndPosition");
        protein_start_position--;
        String seq = hit.getSequence();
        FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(pos, masses, -1, seq, tag);

        if (protein_end_position >= 0) seq = seq.substr(0, protein_end_position);
        if (protein_start_position >= 0) seq = seq.substr(protein_start_position);
        FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(pos_in_truncated, masses_in_truncated, flanking_mass_tol, seq, tag);
        if (!pos_in_truncated.empty())
        {
          for (int i = 0; i < pos.size(); i++)
          {
            if (pos[i] - (protein_start_position < 0 ? 0 : protein_start_position) != pos_in_truncated[0]) continue;
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
      for (size_t i = 0; i <= tag.getLength(); i++) // Fixed signed/unsigned comparison issue
      {
        ss << std::to_string(tag.getScore(i)) << ",";
      }
      ss << "\n";
    }
  }
  fs << ss.str();
}

String FLASHTnTFile::generateProFormaString_(const String& sequence,
                                             int seq_start,
                                             int seq_end,
                                             const std::vector<double>& mod_masses,
                                             const std::vector<int>& mod_starts,
                                             const std::vector<int>& mod_ends,
                                             const std::vector<String>& mod_ids)
{
  if (seq_start < 0) seq_start = 1;
  if (seq_end < 0) seq_end = sequence.length();
  const auto truncated_seq = sequence.substr(seq_start, seq_end - seq_start);
  ProForma proforma(AASequence::fromString(truncated_seq)); // Create ProForma object

  // Loop through modifications and add them to the sequence
  for (size_t i = 0; i < mod_masses.size(); ++i)
  {
    proforma.addModification(mod_starts[i] - seq_start + 1, mod_ends[i] - seq_start + 1, mod_ids[i], mod_masses[i]);
  }

  // Build the ProForma string with range modifications
  String proforma_str = proforma.toProFormaString();

  return proforma_str; // Return the final ProForma string
}


void OpenMS::FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs)
{
  std::stringstream ss;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    String tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists("TagIndices")) indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

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

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma for sequence generation
    String proformaStr = generateProFormaString_(hit.getSequence(), start_in_seq, end_in_seq, mod_masses, mod_starts, mod_ends, mod_ids);
    ss << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass")  << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" <<  hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore")  << "\t" << hit.getMetaValue("PrecursorSNR")  << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
  fs << ss.str();
}

void OpenMS::FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr)
{
  std::stringstream ss;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    if (! hit.metaValueExists("Representative")) continue;
    if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
    String tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists("TagIndices")) indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

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

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma to generate ProForma string
    String proformaStr = generateProFormaString_(hit.getSequence(), start_in_seq, end_in_seq, mod_masses, mod_starts, mod_ends, mod_ids);

    ss << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass") << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" << hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore") << "\t" << hit.getMetaValue("PrecursorSNR") << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
  fs << ss.str();
}

} // namespace OpenMS
// namespace OpenMS