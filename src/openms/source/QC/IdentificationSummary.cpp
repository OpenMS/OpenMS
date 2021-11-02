// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  IdentificationSummary::Result IdentificationSummary::compute(vector<ProteinIdentification>& prot_ids,
                                                               vector<PeptideIdentification>& pep_ids)
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
        if (temp_hits.empty()) continue;
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
    return peptide_spectrum_matches == rhs.peptide_spectrum_matches
          && unique_peptides.count == rhs.unique_peptides.count
          && unique_peptides.fdr_threshold == rhs.unique_peptides.fdr_threshold
          && unique_proteins.count == rhs.unique_proteins.count
          && unique_proteins.fdr_threshold == rhs.unique_proteins.fdr_threshold
          && missed_cleavages_mean == rhs.missed_cleavages_mean
          && protein_hit_scores_mean == rhs.protein_hit_scores_mean
          && peptide_length_mean == rhs.peptide_length_mean;
  }
 

  /// Returns the name of the metric
  const String& IdentificationSummary::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status IdentificationSummary::requires() const
  {
    return QCBase::Status(QCBase::Requires::ID);
  }
}
