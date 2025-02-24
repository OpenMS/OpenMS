// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>


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

  PeptideIdentification::~PeptideIdentification() noexcept = default;

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
    return !std::isnan(rt_);
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
    return !std::isnan(mz_);
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

  void PeptideIdentification::insertHit(PeptideHit&& hit)
  {
    hits_.push_back(std::move(hit));
  }

  void PeptideIdentification::setHits(const std::vector<PeptideHit>& hits)
  {
    hits_ = hits;
  }

  void PeptideIdentification::setHits(std::vector<PeptideHit>&& hits)
  {
    hits_ = std::move(hits);
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

  String PeptideIdentification::getSpectrumReference() const
  {
    return this->getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, "");
  }

  void PeptideIdentification::setSpectrumReference(const String& id)
  {
    this->setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, id);
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
    return this->getMetaValue("experiment_label", "");
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
    return id_.empty()
           && hits_.empty()
           && significance_threshold_ == 0.0
           && score_type_.empty()
           && higher_score_better_ == true
           && base_name_.empty();
  }

  std::vector<PeptideHit> PeptideIdentification::getReferencingHits(const std::vector<PeptideHit>& hits, const std::set<String>& accession)
  {
    std::vector<PeptideHit> filtered;
    for (const PeptideHit& h_it : hits)
    {
      set<String> hit_accessions = h_it.extractProteinAccessionsSet();
      set<String> intersect;
      set_intersection(hit_accessions.begin(), hit_accessions.end(), accession.begin(), accession.end(), std::inserter(intersect, intersect.begin()));
      if (!intersect.empty())
      {
        filtered.push_back(h_it);
      }
    }
    return filtered;
  }

  std::multimap<String, std::pair<Size, Size>>
  PeptideIdentification::buildUIDsFromAllPepIDs(const ConsensusMap &cmap)
  {
    multimap<String, std::pair<Size, Size>> customID_to_cpepID{};

    ProteinIdentification::Mapping mp_c(cmap.getProteinIdentifications());
    //Iterates of the vector of PeptideIdentification to build the UID
    //and the pep_index
    auto lamda = [](const vector<PeptideIdentification> &cpep_ids,
                    const map<String, StringList> &identifier_to_msrunpath,
                    multimap<String, std::pair<Size, Size>> &customID_to_cpepID,
                    const Size &cf_index) {
        Size pep_index = 0;
        for (const PeptideIdentification &cpep_id : cpep_ids)
        {
          std::pair<Size, Size> index = {cf_index, pep_index};
          auto uid = buildUIDFromPepID(cpep_id, identifier_to_msrunpath);
          customID_to_cpepID.insert(make_pair(uid, index));
          ++pep_index;
        }
    };
    //Build the multimap of all UIDs from the ConsensusMap
    //The multimap maps the UID of the PeptideIdentification to all occurrences in the ConsensusMap
    //An occurrence is described as a pair of an index of the ConsensusFeature and the PeptideIdentification
    for (Size i = 0; i < cmap.size(); ++i)
    {
      lamda(cmap[i].getPeptideIdentifications(), mp_c.identifier_to_msrunpath, customID_to_cpepID, i);
    }
    //all unassigned PeptideIdentifications get -1 for the ConsensusFeature index
    lamda(cmap.getUnassignedPeptideIdentifications(), mp_c.identifier_to_msrunpath, customID_to_cpepID, -1);

    return customID_to_cpepID;
  }

  String PeptideIdentification::buildUIDFromPepID(const PeptideIdentification &pep_id,
                                                  const std::map<String, StringList> &fidentifier_to_msrunpath)
  {
    if (!pep_id.metaValueExists("spectrum_reference"))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Spectrum reference missing at PeptideIdentification.");
    }
    String UID;
    const auto &ms_run_path = fidentifier_to_msrunpath.at(pep_id.getIdentifier());
    if (ms_run_path.size() == 1)
    {
      UID = ms_run_path[0] + '|' + pep_id.getSpectrumReference();
    }
    else if (pep_id.metaValueExists("map_index"))
    {
      UID = pep_id.getMetaValue("map_index").toString() + '|' +
            pep_id.getSpectrumReference();
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Multiple files in a run, but no map_index in PeptideIdentification found.");
    }
    return UID;
  }

  
} // namespace OpenMS
