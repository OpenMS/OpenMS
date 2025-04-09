// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/SYSTEM/File.h>

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
    states(BaseFeature::SIZE_OF_ANNOTATIONSTATE, 0) // initialize all with 0
  {
  }

  AnnotationStatistics::AnnotationStatistics(const AnnotationStatistics& rhs) = default;

  AnnotationStatistics& AnnotationStatistics::operator=(const AnnotationStatistics& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }
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
    MetaInfoInterface(),
    RangeManagerContainerType(),
    DocumentIdentifier(),
    ExposedVector<Feature>(),
    UniqueIdInterface(),
    UniqueIdIndexer<FeatureMap>(),
    protein_identifications_(),
    unassigned_peptide_identifications_(),
    data_processing_(),
    id_data_()
  {
  }

  FeatureMap::FeatureMap(const FeatureMap& source) :
    MetaInfoInterface(source),
    RangeManagerContainerType(source),
    DocumentIdentifier(source),
    ExposedVector<Feature>(source),
    UniqueIdInterface(source),
    UniqueIdIndexer<FeatureMap>(source),
    protein_identifications_(source.protein_identifications_),
    unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
    data_processing_(source.data_processing_),
    id_data_() // updated below
  {
    // copy ID data and update references in features:
    IdentificationData::RefTranslator trans = id_data_.merge(source.id_data_);
    for (Feature& feature : *this)
    {
      feature.updateAllIDReferences(trans);
    }
  }

  FeatureMap::FeatureMap(FeatureMap&& source) = default;

  FeatureMap::~FeatureMap() = default;

  FeatureMap& FeatureMap::operator=(const FeatureMap& rhs)  // TODO: cannot be defaulted since OpenMS::IdentificationData is missing operator=
  {
    if (&rhs == this)
    {
      return *this;
    }
    MetaInfoInterface::operator=(rhs);
    RangeManagerType::operator=(rhs);
    DocumentIdentifier::operator=(rhs);
    UniqueIdInterface::operator=(rhs);
    data_ = rhs.data_;
    protein_identifications_ = rhs.protein_identifications_;
    unassigned_peptide_identifications_ = rhs.unassigned_peptide_identifications_;
    data_processing_ = rhs.data_processing_;

    // copy ID data and update references in features:
    id_data_.clear();
    IdentificationData::RefTranslator trans = id_data_.merge(rhs.id_data_);
    for (Feature& feature : *this)
    {
      feature.updateAllIDReferences(trans);
    }

    return *this;
  }

  //FeatureMap& FeatureMap::operator=(FeatureMap&&) = default; // TODO: cannot be defaulted since OpenMS::IdentificationData is missing operator=


  bool FeatureMap::operator==(const FeatureMap& rhs) const
  {
    return data_ == rhs.data_ &&
           MetaInfoInterface::operator==(rhs) &&
           RangeManagerType::operator==(rhs) &&
           DocumentIdentifier::operator==(rhs) &&
           UniqueIdInterface::operator==(rhs) &&
           protein_identifications_ == rhs.protein_identifications_ &&
           unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ &&
           data_processing_ == rhs.data_processing_;
    // @TODO: implement "operator==" for IdentificationData?
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

    if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty())
    {
      OPENMS_LOG_INFO << "DocumentIdentifiers are lost during merge of FeatureMaps\n";
    }
    DocumentIdentifier::operator=(empty_map);

    UniqueIdInterface::operator=(empty_map);

    // merge these:
    protein_identifications_.insert(protein_identifications_.end(), rhs.protein_identifications_.begin(), rhs.protein_identifications_.end());
    unassigned_peptide_identifications_.insert(unassigned_peptide_identifications_.end(), rhs.unassigned_peptide_identifications_.begin(), rhs.unassigned_peptide_identifications_.end());
    data_processing_.insert(data_processing_.end(), rhs.data_processing_.begin(), rhs.data_processing_.end());

    Size n_old_features = size();
    // append features:
    this->insert(this->end(), rhs.begin(), rhs.end());

    // todo: check for double entries
    // features, unassignedpeptides, proteins...

    // merge IDs (new format):
    IdentificationData::RefTranslator trans = id_data_.merge(rhs.id_data_);
    // update ID references of new features:
    for (Size i = n_old_features; i < size(); ++i)
    {
      operator[](i).updateAllIDReferences(trans);
    }

    // consistency
    try
    {
      UniqueIdIndexer<FeatureMap>::updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition&) // assign new UID's for conflicting entries
    {
      Size replaced_uids = UniqueIdIndexer<FeatureMap>::resolveUniqueIdConflicts();
      OPENMS_LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
    }

    return *this;
  }

  void FeatureMap::sortByIntensity(bool reverse)
  {
    if (reverse)
    {
      std::sort(this->begin(), this->end(), [](auto &left, auto &right) {Feature::IntensityLess cmp; return cmp(right, left);});
    }
    else
    {
      std::sort(this->begin(), this->end(), Feature::IntensityLess());
    }
  }

  void FeatureMap::sortByPosition()
  {
    std::sort(this->begin(), this->end(), Feature::PositionLess());
  }

  void FeatureMap::sortByRT()
  {
    std::sort(this->begin(), this->end(), Feature::RTLess());
  }

  void FeatureMap::sortByMZ()
  {
    std::sort(this->begin(), this->end(), Feature::MZLess());
  }

  void FeatureMap::sortByOverallQuality(bool reverse)
  {
    if (reverse)
    {
      std::sort(this->begin(), this->end(), [](auto& left, auto& right) {Feature::OverallQualityLess cmp; return cmp(right, left);});
    }
    else
    {
      std::sort(this->begin(), this->end(), Feature::OverallQualityLess());
    }
  }

  void FeatureMap::updateRanges()
  {
    clearRanges();
    for (const auto& f : *this)
    {
      extendRT(f.getRT());
      extendMZ(f.getMZ());
      extendIntensity(f.getIntensity());
    }

    // enlarge the range by the convex hull points
    for (Size i = 0; i < this->size(); ++i)
    {
      const DBoundingBox<2>& box = this->operator[](i).getConvexHull().getBoundingBox();
      if (!box.isEmpty())
      {
        extendRT(box.minPosition()[Peak2D::RT]);
        extendRT(box.maxPosition()[Peak2D::RT]);
        extendMZ(box.minPosition()[Peak2D::MZ]);
        extendMZ(box.maxPosition()[Peak2D::MZ]);
      }
    }
  }

  void FeatureMap::swapFeaturesOnly(FeatureMap& from)
  {
    data_.swap(from.data_);

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
    id_data_.swap(from.id_data_);
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

  /// set the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
  void FeatureMap::setPrimaryMSRunPath(const StringList& s)
  {
    if (s.empty())
    {
      OPENMS_LOG_WARN << "Setting empty MS runs paths." << std::endl;
      this->setMetaValue("spectra_data", DataValue(s));
      return;
    }

    for (const String& filename : s)
    {
      if (!filename.hasSuffix("mzML") && !filename.hasSuffix("mzml"))
      {
        OPENMS_LOG_WARN << "To ensure tracability of results please prefer mzML files as primary MS run." << std::endl
                        << "Filename: '" << filename << "'" << std::endl;
      }
    }

    this->setMetaValue("spectra_data", DataValue(s));
  }


  void FeatureMap::setPrimaryMSRunPath(const StringList& s, MSExperiment& e)
  {
    StringList ms_path;
    e.getPrimaryMSRunPath(ms_path);
    if (ms_path.size() == 1 && ms_path[0].hasSuffix("mzML") && File::exists(ms_path[0]))
    {
      setPrimaryMSRunPath(ms_path);
    }
    else
    {
      setPrimaryMSRunPath(s);
    }
  }


  /// get the file path to the first MS run
  void FeatureMap::getPrimaryMSRunPath(StringList& toFill) const
  {
    if (this->metaValueExists("spectra_data"))
    {
      toFill = this->getMetaValue("spectra_data");
    }

    if (toFill.empty())
    {
      OPENMS_LOG_WARN << "No MS run annotated in feature map. Setting to 'UNKNOWN' " << std::endl;
      toFill.push_back("UNKNOWN");
    }
  }

  void FeatureMap::clear(bool clear_meta_data)
  {
    data_.clear();

    if (clear_meta_data)
    {
      clearMetaInfo();
      clearRanges();
      this->DocumentIdentifier::operator=(DocumentIdentifier()); // no "clear" method
      clearUniqueId();
      protein_identifications_.clear();
      unassigned_peptide_identifications_.clear();
      data_processing_.clear();
      id_data_.clear();
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


  std::set<IdentificationDataInternal::ObservationMatchRef> FeatureMap::getUnassignedIDMatches() const
  {
    std::set<IdentificationData::ObservationMatchRef> all_matches;
    for (auto it = id_data_.getObservationMatches().begin();
         it != id_data_.getObservationMatches().end(); ++it)
    {
      all_matches.insert(it);
    }
    std::set<IdentificationData::ObservationMatchRef> assigned_matches;
    for (const Feature& feat : *this)
    {
      assigned_matches.insert(feat.getIDMatches().begin(), feat.getIDMatches().end());
      // @TODO: consider subordinate features? - probably not
    }
    std::set<IdentificationData::ObservationMatchRef> result;
    std::set_difference(all_matches.begin(), all_matches.end(),
                        assigned_matches.begin(), assigned_matches.end(),
                        inserter(result, result.end()));
    return result;
  }


  const IdentificationData& FeatureMap::getIdentificationData() const
  {
    return id_data_;
  }


  IdentificationData& FeatureMap::getIdentificationData()
  {
    return id_data_;
  }

}
