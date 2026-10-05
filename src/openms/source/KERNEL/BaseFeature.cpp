// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  const std::string BaseFeature::NamesOfAnnotationState[] =
    {"no ID", "single ID", "multiple IDs (identical)", "multiple IDs (divergent)"};


  BaseFeature::BaseFeature() :
    RichPeak2D(), quality_(0.0), charge_(0), width_(0)
  {
  }

  BaseFeature::BaseFeature(const BaseFeature& rhs, UInt64 map_index):
      RichPeak2D(rhs),
      quality_(rhs.quality_),
      charge_(rhs.charge_),
      width_(rhs.width_),
      peptides_(rhs.peptides_),
      primary_id_(rhs.primary_id_),
      id_matches_(rhs.id_matches_),
      id_queries_(rhs.id_queries_)
  {
    for (auto& pep : this->peptides_)
    {
      pep.setMetaValue("map_index", map_index);
    }
  }

  BaseFeature::BaseFeature(const RichPeak2D& point) :
    RichPeak2D(point), quality_(0.0), charge_(0), width_(0)
  {
  }

  BaseFeature::BaseFeature(const FeatureHandle& fh) :
    RichPeak2D(fh),
    quality_(0.0),
    charge_(fh.getCharge()),
    width_(fh.getWidth()),
    peptides_()
  {
  }

  BaseFeature::BaseFeature(const Peak2D& point) :
    RichPeak2D(point), quality_(0.0), charge_(0), width_(0)
  {
  }

  bool BaseFeature::operator==(const BaseFeature& rhs) const
  {
    return RichPeak2D::operator==(rhs) && (quality_ == rhs.quality_) && (charge_ == rhs.charge_) && (width_ == rhs.width_)
           && (peptides_ == rhs.peptides_) && (primary_id_ == rhs.primary_id_) && (id_matches_ == rhs.id_matches_)
           && (id_queries_ == rhs.id_queries_);
  }

  bool BaseFeature::operator!=(const BaseFeature& rhs) const
  {
    return !operator==(rhs);
  }

  BaseFeature::~BaseFeature() = default;

  BaseFeature::QualityType BaseFeature::getQuality() const
  {
    return quality_;
  }

  void BaseFeature::setQuality(BaseFeature::QualityType quality)
  {
    quality_ = quality;
  }

  BaseFeature::WidthType BaseFeature::getWidth() const
  {
    return width_;
  }

  void BaseFeature::setWidth(BaseFeature::WidthType fwhm)
  {
    // !!! Dirty hack: as long as featureXML doesn't support a width field,
    // we abuse the meta information for this.
    // See also FeatureXMLFile::readFeature_().
    width_ = fwhm;
    setMetaValue("FWHM", fwhm);
  }

  const BaseFeature::ChargeType& BaseFeature::getCharge() const
  {
    return charge_;
  }

  void BaseFeature::setCharge(const BaseFeature::ChargeType& charge)
  {
    charge_ = charge;
  }

  const PeptideIdentificationList& BaseFeature::getPeptideIdentifications()
  const
  {
    return peptides_;
  }

  PeptideIdentificationList& BaseFeature::getPeptideIdentifications()
  {
    return peptides_;
  }

  void BaseFeature::setPeptideIdentifications(
    const PeptideIdentificationList& peptides)
  {
    peptides_ = peptides;
  }

  void BaseFeature::sortPeptideIdentifications()
  {
    // Sort the hits of every identification in a pass of its own. Doing it inside the comparator
    // skips every identification the sort never compares (e.g. the only identification of a
    // feature) and lets the comparator modify its arguments.
    for (PeptideIdentification& pep : peptides_)
    {
      pep.sort();
    }

    // Read the score orientation once, from the first identification that has hits (default:
    // higher is better), and use it for every comparison. Taking it from the left operand makes
    // the comparator asymmetric as soon as two identifications disagree on isHigherScoreBetter(),
    // and the sort requires a strict weak ordering.
    bool higher_score_better = true;
    for (const PeptideIdentification& pep : peptides_)
    {
      if (!pep.getHits().empty())
      {
        higher_score_better = pep.isHigherScoreBetter();
        break;
      }
    }

    // Best first: the identification whose first (= best) hit scores best comes first,
    // identifications without hits go last, and ties keep their relative order (stable sort).
    // Hit presence is the guard, not PeptideIdentification::empty(): an identification with an
    // identifier or score type but no hits is not empty(), yet it has no first hit to read.
    std::stable_sort(peptides_.begin(), peptides_.end(),
                     [higher_score_better](const PeptideIdentification& p1, const PeptideIdentification& p2)
                     {
                       if (p1.getHits().empty())
                       {
                         return false; // no hits: never precedes anything
                       }
                       if (p2.getHits().empty())
                       {
                         return true; // hits precede no hits
                       }
                       const double s1 = p1.getHits()[0].getScore();
                       const double s2 = p2.getHits()[0].getScore();
                       return higher_score_better ? (s1 > s2) : (s1 < s2);
                     });
  }

  BaseFeature::AnnotationState BaseFeature::getAnnotationState() const
  {
    if (id_matches_.empty()) // consider IDs in old format
    {
      if (peptides_.empty())
      {
        return AnnotationState::FEATURE_ID_NONE;
      }
      if (peptides_.size() == 1 && !peptides_[0].getHits().empty())
      {
        return AnnotationState::FEATURE_ID_SINGLE;
      }
      std::set<std::string> seqs;
      for (Size i = 0; i < peptides_.size(); ++i)
      {
        if (!peptides_[i].getHits().empty())
        {
          PeptideIdentification id_tmp = peptides_[i];
          id_tmp.sort();  // look at best hit only - requires sorting
          seqs.insert(id_tmp.getHits()[0].getSequence().toString());
        }
      }
      if (seqs.size() == 1)
      {
        return AnnotationState::FEATURE_ID_MULTIPLE_SAME; // hits have identical seqs
      }
      if (seqs.size() > 1)
      {
        return AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT; // multiple different annotations ... probably bad mapping
      }
      /*else if (seqs.size()==0)*/
      return AnnotationState::FEATURE_ID_NONE;   // very rare case of empty hits
    }
    else // consider IDs in new format
    {
      if (id_matches_.size() == 1)
      {
        return AnnotationState::FEATURE_ID_SINGLE;
      }
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Comparing multiple owning match annotations requires their IdentificationData");
    }
  }

  BaseFeature::AnnotationState BaseFeature::getAnnotationState(const IdentificationData& data) const
  {
    if (id_matches_.empty()) return getAnnotationState();
    std::optional<IdentificationData::MoleculeIdentity> molecule;
    bool divergent = false;
    for (const auto& reference : id_matches_)
    {
      const auto* run = data.findRunByUuid(reference.run_uuid);
      const auto* match = run ? run->findMatch(reference.match) : nullptr;
      if (! match) throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature association refers to a missing match");
      IdentificationData::MoleculeIdentity identity {match->encoding, match->representation};
      if (molecule && *molecule != identity) divergent = true;
      molecule = std::move(identity);
    }
    if (id_matches_.size() == 1) return AnnotationState::FEATURE_ID_SINGLE;
    return divergent ? AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT : AnnotationState::FEATURE_ID_MULTIPLE_SAME;
  }

  bool BaseFeature::hasPrimaryID() const
  {
    return bool(primary_id_);
  }


  const IdentificationData::MoleculeIdentity& BaseFeature::getPrimaryID() const
  {
    if (!primary_id_)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "no primary ID assigned");
    }

    return *primary_id_; // unpack the option
  }


  void BaseFeature::clearPrimaryID()
  {
    primary_id_ = nullopt;
  }


  void BaseFeature::setPrimaryID(const IdentificationData::MoleculeIdentity& id)
  {
    primary_id_ = id;
  }


  const std::set<IdentificationData::MatchReference>& BaseFeature::getIDMatches() const
  {
    return id_matches_;
  }


  std::set<IdentificationData::MatchReference>& BaseFeature::getIDMatches()
  {
    return id_matches_;
  }


  void BaseFeature::addIDMatch(IdentificationData::MatchReference ref)
  {
    id_matches_.insert(ref);
  }


} // namespace OpenMS
