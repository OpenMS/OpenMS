// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <climits>

using namespace std;

namespace OpenMS
{
  IDFilter::IDFilter()
  {
  }

  IDFilter::~IDFilter()
  {
  }


  struct IDFilter::HasMinPeptideLength
  {
    typedef PeptideHit argument_type; // for use as a predicate

    Size length;

    HasMinPeptideLength(Size length):
      length(length)
    {}
    
    bool operator()(const PeptideHit& hit) const
    {
      return hit.getSequence().size() >= length;
    }
  };


  struct IDFilter::HasMinCharge
  {
    typedef PeptideHit argument_type; // for use as a predicate

    Int charge;

    HasMinCharge(Int charge):
      charge(charge)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getCharge() >= charge;
    }
  };


  struct IDFilter::HasLowMZError
  {
    typedef PeptideHit argument_type; // for use as a predicate

    double precursor_mz, tolerance;

    HasLowMZError(double precursor_mz, double tolerance, bool unit_ppm):
      precursor_mz(precursor_mz), tolerance(tolerance)
    {
      if (unit_ppm) tolerance *= precursor_mz / 1.0e6;
    }

    bool operator()(const PeptideHit& hit) const
    {
      Int z = hit.getCharge();
      if (z == 0) z = 1;
      double peptide_mz = (hit.getSequence().getMonoWeight(Residue::Full, z) /
                           double(z));
      return fabs(precursor_mz - peptide_mz) <= tolerance;
    }
  };


  struct IDFilter::HasMatchingModification
  {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& mods;

    HasMatchingModification(const set<String>& mods):
      mods(mods)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      const AASequence& seq = hit.getSequence();
      if (mods.empty()) return seq.isModified();

      for (Size i = 0; i < seq.size(); ++i)
      {
        if (seq.isModified(i))
        {
          String mod_name = seq[i].getModification() + " (" +
            seq[i].getOneLetterCode() + ")";
          if (mods.count(mod_name) > 0) return true;
        }
      }

      // terminal modifications:
      if (seq.hasNTerminalModification())
      {
        String mod_name = seq.getNTerminalModification() + " (N-term)";
        if (mods.count(mod_name) > 0) return true;
      }
      if (seq.hasCTerminalModification())
      {
        String mod_name = seq.getCTerminalModification() + " (C-term)";
        if (mods.count(mod_name) > 0) return true;
      }

      return false;
    }
  };


  struct IDFilter::HasMatchingSequence
  {
    typedef PeptideHit argument_type; // for use as a predicate

    const set<String>& sequences;
    bool ignore_mods;

    HasMatchingSequence(const set<String>& sequences, bool ignore_mods):
      sequences(sequences), ignore_mods(ignore_mods)
    {}

    bool operator()(const PeptideHit& hit) const
    {
      const String& query = (ignore_mods ? 
                             hit.getSequence().toUnmodifiedString() :
                             hit.getSequence().toString());
      return (sequences.count(query) > 0);
    }
  };


  struct IDFilter::HasNoEvidence
  {
    typedef PeptideHit argument_type; // for use as a predicate

    bool operator()(const PeptideHit& hit) const
    {
      return hit.getPeptideEvidences().empty();
    }
  };


  struct IDFilter::HasRTInRange
  {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double rt_min, rt_max;

    HasRTInRange(double rt_min, double rt_max):
      rt_min(rt_min), rt_max(rt_max)
    {}

    bool operator()(const PeptideIdentification& id) const
    {
      double rt = id.getRT();
      return (rt >= rt_min) && (rt <= rt_max);
    }
  };


  struct IDFilter::HasMZInRange
  {
    typedef PeptideIdentification argument_type; // for use as a predicate

    double mz_min, mz_max;

    HasMZInRange(double mz_min, double mz_max):
      mz_min(mz_min), mz_max(mz_max)
    {}

    bool operator()(const PeptideIdentification& id) const
    {
      double mz = id.getMZ();
      return (mz >= mz_min) && (mz <= mz_max);
    }
  };


  void IDFilter::filterIdentificationsUnique(
    const PeptideIdentification& identification,
    PeptideIdentification& filtered_identification)
  {
    // there's no "PeptideHit::operator<" defined, so we can't use a set nor
    // "sort" + "unique" from the standard library
    vector<PeptideHit> hits;
    filtered_identification = identification;
    vector<PeptideHit> temp_hits = identification.getHits();

    for (vector<PeptideHit>::iterator it = temp_hits.begin();
         it != temp_hits.end(); ++it)
    {
      if (find(hits.begin(), hits.end(), *it) == hits.end())
      {
        hits.push_back(*it);
      }
    }
    filtered_identification.setHits(hits);
  }


  void IDFilter::filterIdentificationsByMzError(
    const PeptideIdentification& identification, double mass_error,
    bool unit_ppm, PeptideIdentification& filtered_identification)
  {
    filtered_identification = identification;

    struct HasLowMZError error_filter(filtered_identification.getMZ(),
                                      mass_error, unit_ppm);
    keepMatchingItems(filtered_identification.getHits(), error_filter);

    filtered_identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByBestHits(
    const PeptideIdentification& identification,
    PeptideIdentification& filtered_identification, bool strict)
  {
    filtered_identification = identification;
    vector<PeptideHit>& hits = filtered_identification.getHits();
    if (hits.size() > 1)
    {
      filtered_identification.sort();
      double top_score = hits[0].getScore();
      bool higher_better = filtered_identification.isHigherScoreBetter();
      struct HasGoodScore<PeptideHit> good_score(top_score, higher_better);
      if (strict) // only one best score allowed
      {
        if (good_score(hits[1])) // two (or more) best-scoring hits
        {
          hits.clear();
        }
        else
        {
          hits.resize(1);
        }
      }
      else
      {
        // we could use keepMatchingHits() here, but it would be less efficient
        // (since the hits are already sorted by score):
        for (vector<PeptideHit>::iterator it = ++hits.begin(); it != hits.end();
             ++it)
        {
          if (!good_score(*it))
          {
            hits.erase(it, hits.end());
            break;
          }
        }
      }

      filtered_identification.assignRanks();
    }
  }


  void IDFilter::filterIdentificationsByLength(
    const PeptideIdentification& identification,
    PeptideIdentification& filtered_identification, Size min_length,
    Size max_length)
  {
    filtered_identification = identification;
    if (min_length > 0)
    {
      struct HasMinPeptideLength length_filter(min_length);
      keepMatchingItems(filtered_identification.getHits(), length_filter);
    }
    ++max_length; // the predicate tests for ">=", we need ">"
    if (max_length > min_length)
    {
      struct HasMinPeptideLength length_filter(max_length);
      removeMatchingItems(filtered_identification.getHits(), length_filter);
    }

    filtered_identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByCharge(
    const PeptideIdentification& identification, Int min_charge,
    PeptideIdentification& filtered_identification)
  {
    filtered_identification = identification;

    struct HasMinCharge charge_filter(min_charge);
    keepMatchingItems(filtered_identification.getHits(), charge_filter);
    filtered_identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByVariableModifications(
    const PeptideIdentification& identification, 
    vector<String>& fixed_modifications, 
    PeptideIdentification& filtered_identification)
  {
    filtered_identification = identification;

    set<String> selected_mods;
    if (!fixed_modifications.empty())
    {
      vector<String> all_mods;
      ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
      sort(all_mods.begin(), all_mods.end());
      sort(fixed_modifications.begin(), fixed_modifications.end());
    
      set_difference(all_mods.begin(), all_mods.end(),
                     fixed_modifications.begin(), fixed_modifications.end(),
                     inserter(selected_mods, selected_mods.begin()));
    }

    struct HasMatchingModification mod_filter(selected_mods);
    keepMatchingItems(filtered_identification.getHits(), mod_filter);

    filtered_identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByExclusionPeptides(
    const PeptideIdentification& identification,
    const set<String>& peptides, bool ignore_modifications,
    PeptideIdentification& filtered_identification)
  {
    filtered_identification = identification;

    struct HasMatchingSequence seq_filter(peptides, ignore_modifications);
    removeMatchingItems(filtered_identification.getHits(), seq_filter);

    filtered_identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByPValues(
    PeptideIdentification& identification, const String& metavalue_key,
    double p_value, const String& source_tool)
  {
    double cutoff = 1 - p_value;

    // how many hits are missing the meta value?
    Size n_hits = identification.getHits().size();
    struct HasMetaValue<PeptideHit> present_filter(metavalue_key, DataValue());
    keepMatchingItems(identification.getHits(), present_filter);
    Size n_missing = n_hits - identification.getHits().size();
    
    struct HasMaxMetaValue<PeptideHit> max_filter(metavalue_key, 
                                                  DataValue(cutoff));
    keepMatchingItems(identification.getHits(), max_filter);

    if (n_missing > 0)
    {
      LOG_WARN << "Filtering identifications by p-value did not work for "
               << n_missing << " of " << n_hits
               << " hits. Your data is missing a meta-value ('"
               << metavalue_key << "')";
      if (!source_tool.empty()) LOG_WARN << " added by " << source_tool;
      LOG_WARN << "." << endl;
    }

    identification.assignRanks();
  }


  void IDFilter::filterIdentificationsByRTPValues(
    const PeptideIdentification& identification,
    PeptideIdentification& filtered_identification, double p_value)
  {
    filtered_identification = identification;
    filterIdentificationsByPValues(filtered_identification, 
                                   "predicted_RT_p_value", p_value,
                                   "RTPredict");
  }


  void IDFilter::filterIdentificationsByRTFirstDimPValues(
    const PeptideIdentification& identification,
    PeptideIdentification& filtered_identification, double p_value)
  {
    filtered_identification = identification;
    filterIdentificationsByPValues(filtered_identification, 
                                   "predicted_RT_p_value_first_dim", p_value,
                                   "RTPredict");
  }


  void IDFilter::removeUnreferencedProteinHits(
    const ProteinIdentification& identification,
    const vector<PeptideIdentification>& peptide_identifications,
    ProteinIdentification& filtered_identification)
  {
    filtered_identification = identification;
    const String& run_identifier = filtered_identification.getIdentifier();

    // build set of protein accessions that are referenced by peptides:
    set<String> accessions;
    for (vector<PeptideIdentification>::const_iterator pep_it = 
           peptide_identifications.begin(); pep_it != 
           peptide_identifications.end(); ++pep_it)
    {
      // run ID of protein and peptide identification must match:
      if (pep_it->getIdentifier() == run_identifier)
      {
        // extract protein accessions of each peptide hit:
        for (vector<PeptideHit>::const_iterator hit_it =
               pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
             ++hit_it)
        {
          const set<String>& current_accessions = 
            hit_it->extractProteinAccessions();
          accessions.insert(current_accessions.begin(),
                            current_accessions.end());
        }
      }
    }

    // remove all protein hits not referenced by a peptide:
    struct HasMatchingAccession<ProteinHit> acc_filter(accessions);
    keepMatchingItems(filtered_identification.getHits(), acc_filter);
  }


  // rename: "removeReferencesToMissingProteins"?
  void IDFilter::removeUnreferencedPeptideHits(
    const ProteinIdentification& identification,
    vector<PeptideIdentification>& peptide_identifications,
    bool delete_unreferenced_peptide_hits /* = false */)
  {
    const String& run_identifier = identification.getIdentifier();

    // build set of protein accessions
    set<String> accessions;
    for (vector<ProteinHit>::const_iterator hit_it = 
           identification.getHits().begin(); hit_it !=
           identification.getHits().end(); ++hit_it)
    {
      accessions.insert(hit_it->getAccession());
    }

    struct HasMatchingAccession<PeptideEvidence> acc_filter(accessions);

    // remove peptides which are not referenced
    for (vector<PeptideIdentification>::iterator pep_it = 
           peptide_identifications.begin(); pep_it !=
           peptide_identifications.end(); ++pep_it)
    {
      // run id of protein and peptide identification must match
      if (pep_it->getIdentifier() == run_identifier)
      {
        // check protein accessions of each peptide hit
        for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
             hit_it != pep_it->getHits().end(); ++hit_it)
        {
          // no non-const "PeptideHit::getPeptideEvidences" implemented, so we
          // can't use "keepMatchingItems":
          vector<PeptideEvidence> evidences;
          remove_copy_if(hit_it->getPeptideEvidences().begin(),
                         hit_it->getPeptideEvidences().end(),
                         back_inserter(evidences),
                         not1(acc_filter));
          hit_it->setPeptideEvidences(evidences);
        }

        if (delete_unreferenced_peptide_hits)
        {
          removeMatchingItems(pep_it->getHits(), HasNoEvidence());
        }
      }
    }
  }


  void IDFilter::filterIdentificationsByRT(const vector<PeptideIdentification>& identifications, double min_rt, double max_rt, vector<PeptideIdentification>& filtered_identifications)
  {
    filtered_identifications = identifications;

    struct HasRTInRange rt_filter(min_rt, max_rt);
    keepMatchingItems(filtered_identifications, rt_filter);
  }

  void IDFilter::filterIdentificationsByMZ(const vector<PeptideIdentification>& identifications, double min_mz, double max_mz, vector<PeptideIdentification>& filtered_identifications)
  {
    filtered_identifications = identifications;

    struct HasMZInRange mz_filter(min_mz, max_mz);
    keepMatchingItems(filtered_identifications, mz_filter);
  }

  bool IDFilter::updateProteinGroups(
    const vector<ProteinIdentification::ProteinGroup>& groups,
    const vector<ProteinHit>& hits,
    vector<ProteinIdentification::ProteinGroup>& filtered_groups)
  {
    bool valid = true;

    // we'll do lots of look-ups, so use a suitable data structure:
    set<String> valid_accessions;
    for (vector<ProteinHit>::const_iterator hit_it = hits.begin();
         hit_it != hits.end(); ++hit_it)
    {
      valid_accessions.insert(hit_it->getAccession());
    }

    filtered_groups.clear();
    for (vector<ProteinIdentification::ProteinGroup>::const_iterator group_it =
           groups.begin(); group_it != groups.end(); ++group_it)
    {
      ProteinIdentification::ProteinGroup filtered;
      set_intersection(group_it->accessions.begin(), group_it->accessions.end(),
                       valid_accessions.begin(), valid_accessions.end(),
                       inserter(filtered.accessions,
                                filtered.accessions.begin()));
      if (!filtered.accessions.empty())
      {
        if (filtered.accessions.size() < group_it->accessions.size())
        {
          valid = false; // some proteins removed from group
        }
        filtered.probability = group_it->probability;
        filtered_groups.push_back(filtered);
      }
    }
    
    return valid;
  }

} // namespace OpenMS
