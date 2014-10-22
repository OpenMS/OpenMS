// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>

namespace OpenMS
{
  std::ostream& operator<<(std::ostream& os, const AnnotationStatistics& ann)
  {
    os << "Feature annotation with identifications:" << "\n";
    for (Size i = 0; i < ann.states.size(); ++i)
    {
      os << "    " << BaseFeature::NamesOfAnnotationState[i] << ": " << ann.states[i] << "\n";
    }
    os << std::endl;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const FeatureMap& map)
  {
    os << "# -- DFEATUREMAP BEGIN --" << "\n";
    os << "# POS \tINTENS\tOVALLQ\tCHARGE\tUniqueID" << "\n";
    for (FeatureMap::const_iterator iter = map.begin(); iter != map.end(); ++iter)
    {
      os << iter->getPosition() << '\t'
      << iter->getIntensity() << '\t'
      << iter->getOverallQuality() << '\t'
      << iter->getCharge() << '\t'
      << iter->getUniqueId() << "\n";
    }
    os << "# -- DFEATUREMAP END --" << std::endl;
    return os;
  }

  AnnotationStatistics::AnnotationStatistics() :
    states(BaseFeature::SIZE_OF_ANNOTATIONSTATE, 0)     // initialize all with 0
  {
  }

  AnnotationStatistics::AnnotationStatistics(const AnnotationStatistics& rhs) :
    states(rhs.states)
  {
  }

  AnnotationStatistics& AnnotationStatistics::operator=(const AnnotationStatistics& rhs)
  {
    if (this == &rhs) return *this;

    states = rhs.states;
    return *this;
  }

  bool AnnotationStatistics::operator==(const AnnotationStatistics& rhs) const
  {
    return states == rhs.states;
  }

  AnnotationStatistics& AnnotationStatistics::operator+=(BaseFeature::AnnotationState state)
  {
    ++states[static_cast<Size>(state)];
    return *this;
  }

  FeatureMap::FeatureMap() :
    Base(),
    RangeManagerType(),
    DocumentIdentifier(),
    UniqueIdInterface(),
    UniqueIdIndexer<FeatureMap>(),
    protein_identifications_(),
    unassigned_peptide_identifications_(),
    data_processing_()
  {
  }

  FeatureMap::FeatureMap(const FeatureMap& source) :
    Base(source),
    RangeManagerType(source),
    DocumentIdentifier(source),
    UniqueIdInterface(source),
    UniqueIdIndexer<FeatureMap>(source),
    protein_identifications_(source.protein_identifications_),
    unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
    data_processing_(source.data_processing_)
  {
  }

  FeatureMap::~FeatureMap()
  {
  }

  FeatureMap& FeatureMap::operator=(const FeatureMap& rhs)
  {
    if (&rhs == this) return *this;

    Base::operator=(rhs);
    RangeManagerType::operator=(rhs);
    DocumentIdentifier::operator=(rhs);
    UniqueIdInterface::operator=(rhs);
    protein_identifications_ = rhs.protein_identifications_;
    unassigned_peptide_identifications_ = rhs.unassigned_peptide_identifications_;
    data_processing_ = rhs.data_processing_;

    return *this;
  }

  bool FeatureMap::operator==(const FeatureMap& rhs) const
  {
    return std::operator==(*this, rhs) &&
           RangeManagerType::operator==(rhs) &&
           DocumentIdentifier::operator==(rhs) &&
           UniqueIdInterface::operator==(rhs) &&
           protein_identifications_ == rhs.protein_identifications_ &&
           unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ &&
           data_processing_ == rhs.data_processing_;
  }

  bool FeatureMap::operator!=(const FeatureMap& rhs) const
  {
    return !(operator==(rhs));
  }

  FeatureMap FeatureMap::operator+(const FeatureMap& rhs) const
  {
    FeatureMap tmp(*this);
    tmp += rhs;
    return tmp;
  }

  FeatureMap& FeatureMap::operator+=(const FeatureMap& rhs)
  {
    FeatureMap empty_map;
    // reset these:
    RangeManagerType::operator=(empty_map);

    if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty()) LOG_INFO << "DocumentIdentifiers are lost during merge of FeatureMaps\n";
    DocumentIdentifier::operator=(empty_map);

    UniqueIdInterface::operator=(empty_map);

    // merge these:
    protein_identifications_.insert(protein_identifications_.end(), rhs.protein_identifications_.begin(), rhs.protein_identifications_.end());
    unassigned_peptide_identifications_.insert(unassigned_peptide_identifications_.end(), rhs.unassigned_peptide_identifications_.begin(), rhs.unassigned_peptide_identifications_.end());
    data_processing_.insert(data_processing_.end(), rhs.data_processing_.begin(), rhs.data_processing_.end());

    // append features:
    this->insert(this->end(), rhs.begin(), rhs.end());

    // todo: check for double entries
    // features, unassignedpeptides, proteins...

    // consistency
    try
    {
      UniqueIdIndexer<FeatureMap>::updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition /*&e*/) // assign new UID's for conflicting entries
    {
      Size replaced_uids =  UniqueIdIndexer<FeatureMap>::resolveUniqueIdConflicts();
      LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
    }

    return *this;
  }

  void FeatureMap::sortByIntensity(bool reverse)
  {
    if (reverse)
    {
      std::sort(this->begin(), this->end(), reverseComparator(FeatureType::IntensityLess()));
    }
    else
    {
      std::sort(this->begin(), this->end(), FeatureType::IntensityLess());
    }
  }

  void FeatureMap::sortByPosition()
  {
    std::sort(this->begin(), this->end(), FeatureType::PositionLess());
  }

  void FeatureMap::sortByRT()
  {
    std::sort(this->begin(), this->end(), FeatureType::RTLess());
  }

  void FeatureMap::sortByMZ()
  {
    std::sort(this->begin(), this->end(), FeatureType::MZLess());
  }

  void FeatureMap::sortByOverallQuality(bool reverse)
  {
    if (reverse)
    {
      std::sort(this->begin(), this->end(), reverseComparator(FeatureType::OverallQualityLess()));
    }
    else
    {
      std::sort(this->begin(), this->end(), FeatureType::OverallQualityLess());
    }
  }

  void FeatureMap::updateRanges()
  {
    this->clearRanges();
    updateRanges_(this->begin(), this->end());

    //enlarge the range by the convex hull points
    for (Size i = 0; i < this->size(); ++i)
    {
      DBoundingBox<2> box = this->operator[](i).getConvexHull().getBoundingBox();
      if (!box.isEmpty())
      {
        //update RT
        if (box.minPosition()[Peak2D::RT] < this->pos_range_.minPosition()[Peak2D::RT])
        {
          this->pos_range_.setMinX(box.minPosition()[Peak2D::RT]);
        }
        if (box.maxPosition()[Peak2D::RT] > this->pos_range_.maxPosition()[Peak2D::RT])
        {
          this->pos_range_.setMaxX(box.maxPosition()[Peak2D::RT]);
        }
        //update m/z
        if (box.minPosition()[Peak2D::MZ] < this->pos_range_.minPosition()[Peak2D::MZ])
        {
          this->pos_range_.setMinY(box.minPosition()[Peak2D::MZ]);
        }
        if (box.maxPosition()[Peak2D::MZ] > this->pos_range_.maxPosition()[Peak2D::MZ])
        {
          this->pos_range_.setMaxY(box.maxPosition()[Peak2D::MZ]);
        }
      }
    }
  }

  void FeatureMap::swapFeaturesOnly(FeatureMap& from)
  {
    // TODO used by FeatureFinderAlgorithmPicked -- could it also use regular swap?
    Base::swap(from);

    // swap range information (otherwise its false in both maps)
    FeatureMap tmp;
    tmp.RangeManagerType::operator=(* this);
    this->RangeManagerType::operator=(from);
    from.RangeManagerType::operator=(tmp);
  }

  void FeatureMap::swap(FeatureMap& from)
  {
    // swap features and ranges
    swapFeaturesOnly(from);

    // swap DocumentIdentifier
    DocumentIdentifier::swap(from);

    // swap unique id
    UniqueIdInterface::swap(from);

    // swap unique id index
    UniqueIdIndexer<FeatureMap>::swap(from);

    // swap the remaining members
    protein_identifications_.swap(from.protein_identifications_);
    unassigned_peptide_identifications_.swap(from.unassigned_peptide_identifications_);
    data_processing_.swap(from.data_processing_);
  }

  const std::vector<ProteinIdentification>& FeatureMap::getProteinIdentifications() const
  {
    return protein_identifications_;
  }

  std::vector<ProteinIdentification>& FeatureMap::getProteinIdentifications()
  {
    return protein_identifications_;
  }

  void FeatureMap::setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)
  {
    protein_identifications_ = protein_identifications;
  }

  const std::vector<PeptideIdentification>& FeatureMap::getUnassignedPeptideIdentifications() const
  {
    return unassigned_peptide_identifications_;
  }

  std::vector<PeptideIdentification>& FeatureMap::getUnassignedPeptideIdentifications()
  {
    return unassigned_peptide_identifications_;
  }

  void FeatureMap::setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications)
  {
    unassigned_peptide_identifications_ = unassigned_peptide_identifications;
  }

  const std::vector<DataProcessing>& FeatureMap::getDataProcessing() const
  {
    return data_processing_;
  }

  std::vector<DataProcessing>& FeatureMap::getDataProcessing()
  {
    return data_processing_;
  }

  void FeatureMap::setDataProcessing(const std::vector<DataProcessing>& processing_method)
  {
    data_processing_ = processing_method;
  }

  void FeatureMap::clear(bool clear_meta_data)
  {
    Base::clear();

    if (clear_meta_data)
    {
      clearRanges();
      this->DocumentIdentifier::operator=(DocumentIdentifier()); // no "clear" method
      clearUniqueId();
      protein_identifications_.clear();
      unassigned_peptide_identifications_.clear();
      data_processing_.clear();
    }
  }

  AnnotationStatistics FeatureMap::getAnnotationStatistics() const
  {
    AnnotationStatistics result;
    for (ConstIterator iter = this->begin(); iter != this->end(); ++iter)
    {
      result += iter->getAnnotationState();
    }
    return result;
  }

}
