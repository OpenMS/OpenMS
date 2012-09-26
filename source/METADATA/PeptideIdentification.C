// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  PeptideIdentification::PeptideIdentification() :
    MetaInfoInterface(),
    id_(),
    hits_(),
    significance_threshold_(0.0),
    score_type_(),
    higher_score_better_(true)
  {
  }

  PeptideIdentification::PeptideIdentification(const PeptideIdentification & rhs) :
    MetaInfoInterface(rhs),
    id_(rhs.id_),
    hits_(rhs.hits_),
    significance_threshold_(rhs.significance_threshold_),
    score_type_(rhs.score_type_),
    higher_score_better_(rhs.higher_score_better_)
  {
  }

  PeptideIdentification::~PeptideIdentification()
  {
  }

  PeptideIdentification & PeptideIdentification::operator=(const PeptideIdentification & rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    MetaInfoInterface::operator=(rhs);
    id_ = rhs.id_;
    hits_ = rhs.hits_;
    significance_threshold_ = rhs.significance_threshold_;
    score_type_ = rhs.score_type_;
    higher_score_better_ = rhs.higher_score_better_;

    return *this;
  }

  // Equality operator
  bool PeptideIdentification::operator==(const PeptideIdentification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.getHits()
           && significance_threshold_ == rhs.getSignificanceThreshold()
           && score_type_ == rhs.score_type_
           && higher_score_better_ == rhs.higher_score_better_;
  }

  // Inequality operator
  bool PeptideIdentification::operator!=(const PeptideIdentification & rhs) const
  {
    return !(*this == rhs);
  }

  const std::vector<PeptideHit> & PeptideIdentification::getHits() const
  {
    return hits_;
  }

  void PeptideIdentification::insertHit(const PeptideHit & hit)
  {
    hits_.push_back(hit);
  }

  void PeptideIdentification::setHits(const std::vector<PeptideHit> & hits)
  {
    hits_ = hits;
  }

  DoubleReal PeptideIdentification::getSignificanceThreshold() const
  {
    return significance_threshold_;
  }

  void PeptideIdentification::setSignificanceThreshold(DoubleReal value)
  {
    significance_threshold_ = value;
  }

  String PeptideIdentification::getScoreType() const
  {
    return score_type_;
  }

  void PeptideIdentification::setScoreType(const String & type)
  {
    score_type_ = type;
  }

  bool PeptideIdentification::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void PeptideIdentification::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  const String & PeptideIdentification::getIdentifier() const
  {
    return id_;
  }

  void PeptideIdentification::setIdentifier(const String & id)
  {
    id_ = id;
  }

  void PeptideIdentification::assignRanks()
  {
    if (hits_.empty())
    {
      return;
    }
    UInt rank = 1;
    sort();
    vector<PeptideHit>::iterator lit = hits_.begin();
    Real tmpscore = lit->getScore();
    while (lit != hits_.end())
    {
      lit->setRank(rank);
      ++lit;
      if (lit != hits_.end() && lit->getScore() != tmpscore)
      {
        ++rank;
        tmpscore = lit->getScore();
      }
    }
  }

  void PeptideIdentification::sort()
  {
    if (higher_score_better_)
    {
      std::sort(hits_.begin(), hits_.end(), PeptideHit::ScoreMore());
    }
    else
    {
      std::sort(hits_.begin(), hits_.end(), PeptideHit::ScoreLess());
    }
  }

  bool PeptideIdentification::empty() const
  {
    return id_ == ""
           && hits_.empty()
           && significance_threshold_ == 0.0
           && score_type_ == ""
           && higher_score_better_ == true;
  }

  void PeptideIdentification::getReferencingHits(const String & protein_accession, std::vector<PeptideHit> & peptide_hits) const
  {
    vector<String> accession;

    accession.push_back(protein_accession);
    getReferencingHits(accession, peptide_hits);
  }

  void PeptideIdentification::getReferencingHits(const std::vector<String> & accessions, std::vector<PeptideHit> & peptide_hits) const
  {
    for (Size i = 0; i < hits_.size(); ++i)
    {
      vector<String>::const_iterator it = hits_[i].getProteinAccessions().begin();
      while (it != hits_[i].getProteinAccessions().end())
      {
        if (find(accessions.begin(), accessions.end(), *it) != accessions.end())
        {
          peptide_hits.push_back(hits_[i]);
          it = hits_[i].getProteinAccessions().end();
        }
        else
        {
          ++it;
        }
      }
    }
  }

  void PeptideIdentification::getReferencingHits(const std::vector<ProteinHit> & protein_hits, std::vector<PeptideHit> & peptide_hits) const
  {
    vector<String> accessions;

    for (vector<ProteinHit>::const_iterator it = protein_hits.begin();
         it != protein_hits.end();
         ++it)
    {
      accessions.push_back(it->getAccession());
    }
    getReferencingHits(accessions, peptide_hits);
  }

  void PeptideIdentification::getNonReferencingHits(const String & protein_accession, std::vector<PeptideHit> & peptide_hits) const
  {
    vector<String> accession;

    accession.push_back(protein_accession);
    getNonReferencingHits(accession, peptide_hits);
  }

  void PeptideIdentification::getNonReferencingHits(const std::vector<String> & accessions, std::vector<PeptideHit> & peptide_hits) const
  {
    bool found = false;

    for (Size i = 0; i < hits_.size(); ++i)
    {
      found = false;
      vector<String>::const_iterator it = hits_[i].getProteinAccessions().begin();
      while (it != hits_[i].getProteinAccessions().end())
      {
        if (find(accessions.begin(), accessions.end(), *it) != accessions.end())
        {
          found = true;
        }
        ++it;
      }
      if (!found)
      {
        peptide_hits.push_back(hits_[i]);
      }
    }
  }

  void PeptideIdentification::getNonReferencingHits(const std::vector<ProteinHit> & protein_hits, std::vector<PeptideHit> & peptide_hits) const
  {
    vector<String> accessions;

    for (vector<ProteinHit>::const_iterator it = protein_hits.begin();
         it != protein_hits.end();
         ++it)
    {
      accessions.push_back(it->getAccession());
    }
    getNonReferencingHits(accessions, peptide_hits);
  }

} // namespace OpenMS
