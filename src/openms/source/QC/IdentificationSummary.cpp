// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------


#include <OpenMS/QC/IdentificationSummary.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <iostream>
#include <set>
using namespace std;

namespace OpenMS
{

  IdentificationSummary::Result IdentificationSummary::compute(vector<ProteinIdentification>& prot_ids, vector<PeptideIdentification>& pep_ids)
  {
    IdentificationSummary::Result result;
    set<String> peptides;
    set<String> proteins;

    // PSMs and collect unique peptides in set
    for (const auto& pep_id : pep_ids)
    {
      if (!pep_id.empty())
      {
        result.peptide_spectrum_matches += 1;
        const auto& temp_hits = pep_id.getHits();
        if (temp_hits.empty())
          continue;
        peptides.insert(temp_hits[0].getSequence().toUnmodifiedString());
      }
    }
    // get sum of all peptide length for mean calculation
    UInt peptide_length_sum = 0;
    for (const auto& pep : peptides)
    {
      peptide_length_sum += pep.size();
    }
    result.peptide_length_mean = (double)peptide_length_sum / peptides.size();
    // get missed cleavages
    UInt missed_cleavages = 0;
    UInt pep_count = 0;
    MissedCleavages mc;
    mc.compute(prot_ids, pep_ids);
    for (const auto& m : mc.getResults())
    {
      for (const auto& [key, val] : m)
      {
        missed_cleavages += key * val;
        pep_count += val;
      }
    }
    result.missed_cleavages_mean = (double)missed_cleavages / pep_count;
    // collect unique proteins in sets and scores mean
    double protein_hit_scores_sum = 0;
    UInt protein_hit_count = 0;
    for (const auto& prot_id : prot_ids)
    {
      const auto& temp_hits = prot_id.getHits();
      protein_hit_count += temp_hits.size();
      for (const auto& temp_hit : temp_hits)
      {
        proteins.insert(temp_hit.getAccession());
        protein_hit_scores_sum += temp_hit.getScore();
      }
    }
    result.protein_hit_scores_mean = protein_hit_scores_sum / protein_hit_count;
    // unique peptides and proteins with their significance threshold (always the same in idXML file)
    // get significance threshold if score type is FDR, else -1
    result.unique_peptides.count = peptides.size();
    result.unique_proteins.count = proteins.size();
    if (pep_ids.front().getScoreType() == "FDR")
    {
      result.unique_peptides.fdr_threshold = pep_ids.front().getSignificanceThreshold();
    }
    if (prot_ids.front().getScoreType() == "FDR")
    {
      result.unique_proteins.fdr_threshold = prot_ids.front().getSignificanceThreshold();
    }
    return result;
  }

  bool IdentificationSummary::Result::operator==(const Result& rhs) const
  {
    return peptide_spectrum_matches == rhs.peptide_spectrum_matches && unique_peptides.count == rhs.unique_peptides.count && unique_peptides.fdr_threshold == rhs.unique_peptides.fdr_threshold &&
           unique_proteins.count == rhs.unique_proteins.count && unique_proteins.fdr_threshold == rhs.unique_proteins.fdr_threshold && missed_cleavages_mean == rhs.missed_cleavages_mean &&
           protein_hit_scores_mean == rhs.protein_hit_scores_mean && peptide_length_mean == rhs.peptide_length_mean;
  }


  /// Returns the name of the metric
  const String& IdentificationSummary::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status IdentificationSummary::requirements() const
  {
    return QCBase::Status(QCBase::Requires::ID);
  }
} // namespace OpenMS
