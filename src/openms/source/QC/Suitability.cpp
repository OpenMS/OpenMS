// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/QC/Suitability.h>

using namespace std;

namespace OpenMS
{
  double Suitability::computeSpectraQuality(const MSExperiment& exp, const std::vector<PeptideIdentification>& pep_ids)
  {
    num_ms2 = 0;
    for (auto const& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 2)
      {
        ++num_ms2;
      }
    }

    if (num_ms2 == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS2 spectra found");
    }

    num_novo_seqs = 0;
    set<AASequence> unique_novo;

    for (const auto& pep_id : pep_ids)
    {
      if (pep_id.getHits().empty()) continue;
      ++num_novo_seqs;
      unique_novo.insert(pep_id.getHits()[0].getSequence());
    }

    num_unique_novo_seqs = unique_novo.size();
    spectral_quality = double(num_ms2) / num_novo_seqs;

    return spectral_quality;
  }
  /*
  double Suitability::computeSuitability(const vector<PeptideIdentification>& pepIDs, double FDR = 0.01, double novo_fract = 1, bool no_re_rank = false)
  {
    if (!no_re_rank)
    {
      cut_off = getDecoyCutOff_(pepIDs, novo_fract);
      if (cut_off == DBL_MAX)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not compute decoy cut off. Re-ranking impossible. If you want to ignore this, disable re-ranking.");
      }
    }

    num_top_db = 0;
    num_top_novo = 0;
    num_re_ranked = 0;
    num_interest = 0;

    for (const PeptideIdentification& pep_id : pepIDs)
    {
      const vector<PeptideHit>& hits = pep_id.getHits();
      bool q_value_score = (pep_id.getScoreType() == "q-value");

      if (hits.empty()) continue;

      // sort hits by q-value
      if (q_value_score)
      {
        sort(hits.begin(), hits.end(),
          [](const PeptideHit& a, const PeptideHit& b)
          {
            return a.getScore() < b.getScore();
          });
      }
      else
      {
        if (!hits[0].metaValueExists("q-value"))
        {
          throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
        }

        sort(hits.begin(), hits.end(),
          [](const PeptideHit& a, const PeptideHit& b)
          {
            return float(a.getMetaValue("q-value")) < float(b.getMetaValue("q-value"));
          });
      }


      const PeptideHit& top_hit = hits[0];

      // skip if the top hit is a decoy hit
      if (!top_hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }
      if (top_hit.getMetaValue("target_decoy") == "decoy") continue;

      // skip if top hit is out ouf FDR
      if (scoreHigherThanFDR_(top_hit, FDR, q_value_score)) continue;

      // check if top hit is found in de novo protein
      if (!isNovoHit_(top_hit)) // top hit is db hit
      {
        ++num_top_db;
        continue;
      }

      // find the second target hit, skip all decoy or novo hits inbetween
      const PeptideHit* second_hit = nullptr;
      String target = "target";
      for (UInt i = 1; i < hits.size(); ++i)
      {
        // check for FDR
        if (scoreHigherThanFDR_(hits[i], FDR, q_value_score)) break;

        if (target.find(String(hits[i].getMetaValue("target_decoy"), 0)) == 0) // also check for "target+decoy" value
        {
          // check if hit is novo hit
          if (isNovoHit_(hits[i])) continue;

          second_hit = &hits[i];
          break;
        }
      }
      if (second_hit == nullptr) // no second target hit with given FDR found
      {
        ++num_top_novo;
        continue;
      }

      // second hit is db hit
      ++num_interest;

      // check for re-ranking
      if (no_re_rank)
      {
        ++num_top_novo;
        continue;
      }

      // check for xcorr score
      if (!top_hit.metaValueExists("MS:1002252") || !second_hit->metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      double top_xscore_mw = double(top_hit.getMetaValue("MS:1002252")) / top_hit.getSequence().getMonoWeight();
      double second_xscore_mw = double(second_hit->getMetaValue("MS:1002252")) / second_hit->getSequence().getMonoWeight();
      if (top_xscore_mw - second_xscore_mw <= cut_off)
      {
        ++num_top_db;
        ++num_re_ranked;
      }
      else
      {
        ++num_top_novo;
      }
    }

    suitability = double(num_top_db) / (num_top_db + num_top_novo);

    return suitability;
  }

  map<String, double> Suitability::getData()
  {
    map<String, double> data_map;

    data_map["#novor_seqs"] = double(num_novo_seqs);
    data_map["#MS2_spectra"] = double(num_ms2);
    data_map["#unique_novor_seqs"] = double(num_unique_novo_seqs);
    data_map["spectral_quality"] = spectral_quality;

    data_map["#top_novo_hits"] = double(num_top_novo);
    data_map["#top_db_hits"] = double(num_top_db);
    data_map["#re_ranked"] = double(num_re_ranked);
    data_map["#possible_re_ranks"] = double(num_interest);
    data_map["cut_off"] = cut_off;
    data_map["suitability"] = suitability;

    return data_map;
  }

  double Suitability::getDecoyDiff_(const PeptideIdentification& pep_id)
  {
    double diff = DBL_MAX;

    // get the score of the first two decoy hits
    double decoy_1 = DBL_MAX;
    double decoy_2 = DBL_MAX;
    UInt curr_hit = 0;

    for (const auto& hit : pep_id.getHits())
    {
      if (curr_hit > 10) break;
      ++curr_hit;

      if (!hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }

      if (pep_id.getScoreType() != "q-value")
      {
        if (!hit.metaValueExists("q-value"))
        {
          throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
        }
      }

      if (!hit.metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      if (decoy_1 == DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = hit.getMetaValue("MS:1002252");
        continue;
      }
      if (decoy_1 < DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = hit.getMetaValue("MS:1002252");
        break;
      }
    }

    if (decoy_2 < DBL_MAX) // if there are two decoy hits
    {
      diff = abs(decoy_1 - decoy_2) / pep_id.getHits()[0].getSequence().getMonoWeight(); // normalized by mw
    }

    // if there aren't two decoy hits DBL_MAX is returned
    return diff;
  }

  double Suitability::getDecoyCutOff_(const vector<PeptideIdentification>& pep_ids, double novor_fract)
  {
    // get all decoy diffs of peptide ids with at least two decoy hits
    vector<double> diffs;
    for (const auto& pep_id : pep_ids)
    {
      double diff = getDecoyDiff_(pep_id);
      if (diff < DBL_MAX)
      {
        diffs.push_back(diff);
      }
    }

    if (double(diffs.size()) / pep_ids.size() < 0.2)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'force_no_re_rank' flag to still compute a suitability score."));
    }

    // sort the diffs decreasing and get the (1-novo_fract)*N one
    auto sort_end = diffs.begin() + (1 - novor_fract) * diffs.size();

    partial_sort(diffs.begin(), sort_end + 1, diffs.end(), greater<double>());

    return *sort_end;
  }

  bool Suitability::isNovoHit_(const PeptideHit& hit)
  {
    const set<String> accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
      {
        return false;
      }
    }
    return true;
  }

  bool Suitability::scoreHigherThanFDR_(const PeptideHit& hit, double FDR, bool q_value_score)
  {
    if (q_value_score) // score type is q-value
    {
      if (hit.getScore() > FDR) return true;
      return false;
    }

    if (hit.metaValueExists("q-value")) // look for q-value at metavalues
    {
      if (float(hit.getMetaValue("q-value")) > FDR) return true;
      return false;
    }

    // no q-value found
    throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
  }*/
}
