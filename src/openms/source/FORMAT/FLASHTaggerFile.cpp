// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Tagger output *.tsv file format
     @ingroup FileIO
  **/
  void FLASHTaggerFile::writeTagHeader(std::fstream& fs)
  {
    fs << "TagIndex\tProteinIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPos\tLength\tDeNovoScore\tMasses\tMassScores\t"
          "Scans\n";
  }

  /// write header line for topFD feature file
  void FLASHTaggerFile::writeProteinHeader(std::fstream& fs)
  {
    fs << "ProteinIndex\tProteinAccession\tProteinDescription\tProteinSequence\tMatchedAminoAcidCount\tCoverage(%)"
          "\tProteinScore\tProteinQvalue\tTagCount\tTagIndices\n";
  }

  /// write the features in regular file output
  void FLASHTaggerFile::writeTags(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    for (int n = 0; n < tagger.getProteinHits().size(); n++)
    {
      for (const auto& tag : tagger.getTags())
      {
        auto hits = tagger.getProteinHits(tag);
        if (n < tagger.getProteinHits().size())
        {
          bool found = false;
          for (const auto& hit : hits)
          {
            if (n == tagger.getProteinIndex(hit)) found = true;
          }
          if (! found) continue;
        }

        if (n == tagger.getProteinHits().size() && ! hits.empty()) continue;

        String acc = "";
        String description = "";
        String hitindices = "";
        String positions = "";
        for (const auto& hit : hits)
        {
          if (! acc.empty()) acc += ";";
          if (! description.empty()) description += ";";
          if (! hitindices.empty()) hitindices += ";";
          if (! positions.empty()) positions += ";";
          acc += hit.getAccession();
          description += hit.getDescription();
          hitindices += std::to_string(tagger.getProteinIndex(hit));

          auto seqposition = tagger.getMatchedPositions(hit, tag);
          if (seqposition.size() != 0) { positions += std::to_string(seqposition[0]); }
        }

        fs << tagger.getTagIndex(tag) << "\t" << hitindices << "\t" << acc << "\t" << description << "\t" << tag.getSequence() << "\t"
           << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions << "\t" << tag.getLength() << "\t"
           << tag.getScore() << "\t";

        for (const auto& mz : tag.getMzs())
        {
          fs << std::to_string(mz) << ",";
        }
        fs << "\t";
        for (auto i = 0; i < tag.getLength(); i++)
        {
          fs << std::to_string(tag.getScore(i)) << ",";
        }
        fs << "\t";
        for (auto i = 0; i < tag.getLength(); i++)
        {
          fs << std::to_string(tag.getScan(i)) << ",";
        }
        fs << "\n";
      }
    }
  }


  void FLASHTaggerFile::writeProteins(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    for (const auto& hit : tagger.getProteinHits())
    {
      String tagindices = "";
      int cntr = 0;
      for (const auto& tag : tagger.getTags(hit))
      {
        if (! tagindices.empty()) tagindices += ";";
        tagindices += std::to_string(tagger.getTagIndex(tag));
        cntr++;
      }

      fs << tagger.getProteinIndex(hit) << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getSequence() << "\t"
         << hit.getMetaValue("MatchedAA") << "\t" << 100.0 * hit.getCoverage() << "\t" << hit.getScore() << "\t"
         << std::to_string((double)hit.getMetaValue("qvalue")) << "\t" << cntr << "\t" << tagindices << "\n";
    }
  }
} // namespace OpenMS