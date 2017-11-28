// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <sstream>
#include <iostream>

using namespace std;

namespace OpenMS
{

  PeptideIdentification::PeptideIdentification() :
    MetaInfoInterface(),
    id_(),
    hits_(),
    significance_threshold_(0.0),
    score_type_(),
    higher_score_better_(true),
    base_name_(),
    mz_(std::numeric_limits<double>::quiet_NaN()),
    rt_(std::numeric_limits<double>::quiet_NaN())
  {
  }

  PeptideIdentification::PeptideIdentification(const PeptideIdentification& rhs) :
    MetaInfoInterface(rhs),
    id_(rhs.id_),
    hits_(rhs.hits_),
    significance_threshold_(rhs.significance_threshold_),
    score_type_(rhs.score_type_),
    higher_score_better_(rhs.higher_score_better_),
    base_name_(rhs.base_name_),
    mz_(rhs.mz_),
    rt_(rhs.rt_)
  {
    setExperimentLabel( rhs.getExperimentLabel() );
  }

  PeptideIdentification::~PeptideIdentification()
  {
  }

  PeptideIdentification& PeptideIdentification::operator=(const PeptideIdentification& rhs)
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
    setExperimentLabel( rhs.getExperimentLabel() );
    base_name_ = rhs.base_name_;
    mz_ = rhs.mz_;
    rt_ = rhs.rt_;

    return *this;
  }

  // Equality operator
  bool PeptideIdentification::operator==(const PeptideIdentification& rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.hits_
           && significance_threshold_ == rhs.getSignificanceThreshold()
           && score_type_ == rhs.score_type_
           && higher_score_better_ == rhs.higher_score_better_
           && getExperimentLabel() == rhs.getExperimentLabel()
           && base_name_ == rhs.base_name_
           && (mz_ == rhs.mz_ || (!this->hasMZ() && !rhs.hasMZ())) // might be NaN, so comparing == will always be false
           && (rt_ == rhs.rt_ || (!this->hasRT() && !rhs.hasRT()));// might be NaN, so comparing == will always be false
  }

  // Inequality operator
  bool PeptideIdentification::operator!=(const PeptideIdentification& rhs) const
  {
    return !(*this == rhs);
  }

  double PeptideIdentification::getRT() const
  {
    return rt_;
  }

  void PeptideIdentification::setRT(double rt)
  {
    rt_ = rt;
  }

  bool PeptideIdentification::hasRT() const
  {
    return !boost::math::isnan(rt_);
  }

  double PeptideIdentification::getMZ() const
  {
    return mz_;
  }

  void PeptideIdentification::setMZ(double mz)
  {
    mz_ = mz;
  }

  bool PeptideIdentification::hasMZ() const
  {
    return !boost::math::isnan(mz_);
  }

  const std::vector<PeptideHit>& PeptideIdentification::getHits() const
  {
    return hits_;
  }

  std::vector<PeptideHit>& PeptideIdentification::getHits()
  {
    return hits_;
  }

  void PeptideIdentification::insertHit(const PeptideHit& hit)
  {
    hits_.push_back(hit);
  }

  void PeptideIdentification::setHits(const std::vector<PeptideHit>& hits)
  {
    hits_ = hits;
  }

  double PeptideIdentification::getSignificanceThreshold() const
  {
    return significance_threshold_;
  }

  void PeptideIdentification::setSignificanceThreshold(double value)
  {
    significance_threshold_ = value;
  }

  const String& PeptideIdentification::getScoreType() const
  {
    return score_type_;
  }

  void PeptideIdentification::setScoreType(const String& type)
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

  const String& PeptideIdentification::getIdentifier() const
  {
    return id_;
  }

  void PeptideIdentification::setIdentifier(const String& id)
  {
    id_ = id;
  }

  const String& PeptideIdentification::getBaseName() const
  {
    return base_name_;
  }

  void PeptideIdentification::setBaseName(const String& base_name)
  {
    base_name_ = base_name;
  }

  const String PeptideIdentification::getExperimentLabel() const
  {
    // implement as meta value in order to reduce bloat of PeptideIdentification object
    //  -> this is mostly used for pepxml at the moment which allows each peptide id to belong to a different experiment
    if (metaValueExists("experiment_label"))
    {
      return getMetaValue("experiment_label").toString();
    }
    else
    {
      return "";
    }
  }

  void PeptideIdentification::setExperimentLabel(const String& label)
  {
    // do not store empty label (default value)
    if (!label.empty())
    {
      setMetaValue("experiment_label", label);
    }
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
    double last_score = lit->getScore();
    while (lit != hits_.end())
    {
      if ((double)lit->getScore() != last_score)
      {
        ++rank;
        last_score = lit->getScore();
      }
      lit->setRank(rank);
      ++lit;
    }
  }

  void PeptideIdentification::sort()
  {
    if (higher_score_better_)
    {
      std::stable_sort(hits_.begin(), hits_.end(), PeptideHit::ScoreMore());
    }
    else
    {
      std::stable_sort(hits_.begin(), hits_.end(), PeptideHit::ScoreLess());
    }
  }

  void PeptideIdentification::sortByRank()
  {
    std::sort(hits_.begin(), hits_.end(), PeptideHit::RankLess());
  }

  bool PeptideIdentification::empty() const
  {
    return id_ == ""
           && hits_.empty()
           && significance_threshold_ == 0.0
           && score_type_ == ""
           && higher_score_better_ == true
           && base_name_ == "";
  }

  std::vector<PeptideHit> PeptideIdentification::getReferencingHits(const std::vector<PeptideHit>& hits, const std::set<String>& accession)
  {
    std::vector<PeptideHit> filtered;
    for (std::vector<PeptideHit>::const_iterator h_it = hits.begin(); h_it != hits.end(); ++h_it)
    {
      set<String> hit_accessions = h_it->extractProteinAccessionsSet();
      set<String> intersect;
      set_intersection(hit_accessions.begin(), hit_accessions.end(), accession.begin(), accession.end(), std::inserter(intersect, intersect.begin()));
      if (!intersect.empty())
      {
        filtered.push_back(*h_it);
      }
    }
    return filtered;
  }
  
} // namespace OpenMS
