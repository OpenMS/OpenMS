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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

using namespace std;

namespace OpenMS
{
  const std::string BaseFeature::NamesOfAnnotationState[] =
    {"no ID", "single ID", "multiple IDs (identical)", "multiple IDs (divergent)"};


  BaseFeature::BaseFeature() :
    RichPeak2D(), quality_(0.0), charge_(0), width_(0)
  {
  }

  BaseFeature::BaseFeature(const BaseFeature& rhs, UInt64 map_index) :
      RichPeak2D(rhs), quality_(rhs.quality_), charge_(rhs.charge_), width_(rhs.width_),
      peptides_(rhs.peptides_), primary_id_(rhs.primary_id_), id_matches_(rhs.id_matches_)
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
    return RichPeak2D::operator==(rhs)
           && (quality_ == rhs.quality_)
           && (charge_ == rhs.charge_)
           && (width_ == rhs.width_)
           && (peptides_ == rhs.peptides_)
           && (primary_id_ == rhs.primary_id_)
           && (id_matches_ == rhs.id_matches_);
  }

  bool BaseFeature::operator!=(const BaseFeature& rhs) const
  {
    return !operator==(rhs);
  }

  BaseFeature::~BaseFeature()
  {
  }

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

  const vector<PeptideIdentification>& BaseFeature::getPeptideIdentifications()
  const
  {
    return peptides_;
  }

  vector<PeptideIdentification>& BaseFeature::getPeptideIdentifications()
  {
    return peptides_;
  }

  void BaseFeature::setPeptideIdentifications(
    const vector<PeptideIdentification>& peptides)
  {
    peptides_ = peptides;
  }

  void BaseFeature::sortPeptideIdentifications()
  {
    std::sort(peptides_.rbegin(),peptides_.rend(),
              [](PeptideIdentification& p1, PeptideIdentification& p2)
              {p1.sort();p2.sort();
              if (p1.empty())
              {
                return true;
              }
              if (p2.empty())
              {
                return false;
              }
              if (p1.isHigherScoreBetter())
              {
                return p1.getHits()[0].getScore() < p2.getHits()[0].getScore();
              }
              else
              {
                return p1.getHits()[0].getScore() > p2.getHits()[0].getScore();
              }});
  }

  BaseFeature::AnnotationState BaseFeature::getAnnotationState() const
  {
    if (id_matches_.empty()) // consider IDs in old format
    {
      if (peptides_.empty())
      {
        return FEATURE_ID_NONE;
      }
      if (peptides_.size() == 1 && !peptides_[0].getHits().empty())
      {
        return FEATURE_ID_SINGLE;
      }
      std::set<String> seqs;
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
        return FEATURE_ID_MULTIPLE_SAME; // hits have identical seqs
      }
      if (seqs.size() > 1)
      {
        return FEATURE_ID_MULTIPLE_DIVERGENT; // multiple different annotations ... probably bad mapping
      }
      /*else if (seqs.size()==0)*/
      return FEATURE_ID_NONE;   // very rare case of empty hits
    }
    else // consider IDs in new format
    {
      if (id_matches_.size() == 1)
      {
        return FEATURE_ID_SINGLE;
      }
      // if there are multiple IDs, check if all are equal (to the first):
      auto it = id_matches_.begin();
      IdentificationData::IdentifiedMolecule molecule = (*it)->identified_molecule_var;
      for (++it; it != id_matches_.end(); ++it)
      {
        if ((*it)->identified_molecule_var != molecule)
        {
          return FEATURE_ID_MULTIPLE_DIVERGENT;
        }
      }
      return FEATURE_ID_MULTIPLE_SAME;
    }
  }


  bool BaseFeature::hasPrimaryID() const
  {
    return bool(primary_id_);
  }


  const IdentificationData::IdentifiedMolecule& BaseFeature::getPrimaryID() const
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


  void BaseFeature::setPrimaryID(const IdentificationData::IdentifiedMolecule& id)
  {
    primary_id_ = id;
  }


  const std::set<IdentificationData::ObservationMatchRef>& BaseFeature::getIDMatches() const
  {
    return id_matches_;
  }


  std::set<IdentificationData::ObservationMatchRef>& BaseFeature::getIDMatches()
  {
    return id_matches_;
  }


  void BaseFeature::addIDMatch(IdentificationData::ObservationMatchRef ref)
  {
    id_matches_.insert(ref);
  }

  void BaseFeature::updateIDReferences(const IdentificationData::RefTranslator& trans)
  {
    if (primary_id_ != nullopt) // is feature annotated with a "primary ID"?
    {
      primary_id_ = trans.translate(*primary_id_);
    }
    set<IdentificationData::ObservationMatchRef> matches; // refs. to e.g. PSMs
    matches.swap(id_matches_);
    for (const auto& item : matches)
    {
      id_matches_.insert(trans.translate(item));
    }
  }

} // namespace OpenMS
