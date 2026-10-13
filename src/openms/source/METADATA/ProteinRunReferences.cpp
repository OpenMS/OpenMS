// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinRunReferences.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    using RunIdentifiers = std::unordered_set<std::string>;

    RunIdentifiers runIdentifiers(const std::vector<ProteinIdentification>& runs)
    {
      RunIdentifiers identifiers;
      identifiers.reserve(runs.size());
      for (const auto& run : runs)
      {
        identifiers.insert(run.getIdentifier());
      }
      return identifiers;
    }

    void checkPeptides(const RunIdentifiers& runs, const PeptideIdentificationList& peptides)
    {
      for (const auto& peptide : peptides)
      {
        if (!runs.contains(peptide.getIdentifier()))
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            ProteinRunReferences::missingRunMessage(peptide.getIdentifier()));
        }
      }
    }

    void checkFeature(const RunIdentifiers& runs, const Feature& feature)
    {
      checkPeptides(runs, feature.getPeptideIdentifications());
      for (const auto& subordinate : feature.getSubordinates())
      {
        checkFeature(runs, subordinate);
      }
    }

    /// the protein accessions of each run, by run identifier
    using RunAccessions = std::unordered_map<std::string, std::unordered_set<std::string>>;

    RunAccessions runAccessions(const std::vector<ProteinIdentification>& runs)
    {
      RunAccessions accessions;
      for (const auto& run : runs)
      {
        auto& run_accessions = accessions[run.getIdentifier()];
        for (const auto& hit : run.getHits())
        {
          run_accessions.insert(hit.getAccession());
        }
      }
      return accessions;
    }

    void checkAccessions(const RunAccessions& runs, const PeptideIdentificationList& peptides)
    {
      for (const auto& peptide : peptides)
      {
        const auto run = runs.find(peptide.getIdentifier());
        if (run == runs.end())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            ProteinRunReferences::missingRunMessage(peptide.getIdentifier()));
        }
        for (const auto& hit : peptide.getHits())
        {
          for (const auto& evidence : hit.getPeptideEvidences())
          {
            const std::string& accession = evidence.getProteinAccession();
            if (accession.empty()) continue; // not written
            if (!run->second.contains(accession))
            {
              throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "No accession " + accession + " found in run '" + peptide.getIdentifier() + "' for PSM " +
                hit.getSequence().toString() + "_" + std::to_string(hit.getCharge()) +
                ". Every protein of a peptide evidence needs to be a protein hit of the peptide identification's run.");
            }
          }
        }
      }
    }

    void checkFeatureAccessions(const RunAccessions& runs, const Feature& feature)
    {
      checkAccessions(runs, feature.getPeptideIdentifications());
      for (const auto& subordinate : feature.getSubordinates())
      {
        checkFeatureAccessions(runs, subordinate);
      }
    }
  } // namespace

  std::string ProteinRunReferences::missingRunMessage(const std::string& identifier)
  {
    return "Peptide identification has no matching protein run: '" + identifier +
           "'. Every peptide identification needs the protein identification run (search run) with its identifier, which may have no protein hits.";
  }

  void ProteinRunReferences::check(const std::vector<ProteinIdentification>& runs, const std::string& identifier)
  {
    for (const auto& run : runs)
    {
      if (run.getIdentifier() == identifier) return;
    }
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, missingRunMessage(identifier));
  }

  void ProteinRunReferences::check(const std::vector<ProteinIdentification>& runs, const PeptideIdentificationList& peptides)
  {
    if (peptides.empty()) return;
    checkPeptides(runIdentifiers(runs), peptides);
  }

  void ProteinRunReferences::check(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides)
  {
    if (peptides.empty()) return;
    const RunIdentifiers identifiers = runIdentifiers(runs);
    for (const PeptideIdentification* peptide : peptides)
    {
      if (!identifiers.contains(peptide->getIdentifier()))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, missingRunMessage(peptide->getIdentifier()));
      }
    }
  }

  void ProteinRunReferences::check(const FeatureMap& map)
  {
    const RunIdentifiers runs = runIdentifiers(map.getProteinIdentifications());
    for (const auto& feature : map)
    {
      checkFeature(runs, feature);
    }
    checkPeptides(runs, map.getUnassignedPeptideIdentifications());
  }

  void ProteinRunReferences::check(const ConsensusMap& map)
  {
    const RunIdentifiers runs = runIdentifiers(map.getProteinIdentifications());
    for (const auto& feature : map)
    {
      checkPeptides(runs, feature.getPeptideIdentifications());
    }
    checkPeptides(runs, map.getUnassignedPeptideIdentifications());
  }

  void ProteinRunReferences::checkProteinAccessions(const FeatureMap& map)
  {
    const RunAccessions runs = runAccessions(map.getProteinIdentifications());
    for (const auto& feature : map)
    {
      checkFeatureAccessions(runs, feature);
    }
    checkAccessions(runs, map.getUnassignedPeptideIdentifications());
  }

  void ProteinRunReferences::checkProteinAccessions(const ConsensusMap& map)
  {
    const RunAccessions runs = runAccessions(map.getProteinIdentifications());
    for (const auto& feature : map)
    {
      checkAccessions(runs, feature.getPeptideIdentifications());
    }
    checkAccessions(runs, map.getUnassignedPeptideIdentifications());
  }
} // namespace OpenMS
