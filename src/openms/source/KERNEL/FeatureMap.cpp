// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

  namespace
  {
    /// Replace the links of @p feature to erased matches by links to their identifications (see eraseMatches())
    void relinkErased(BaseFeature& feature, const IdentificationData& data,
                      const std::map<IdentificationData::MatchReference, IdentificationData::QueryReference>& erased)
    {
      std::set<IdentificationData::QueryReference> orphaned;
      auto& matches = feature.getIDMatches();
      for (auto it = matches.begin(); it != matches.end();)
      {
        const auto found = erased.find(*it);
        if (found == erased.end())
        {
          ++it;
          continue;
        }
        orphaned.insert(found->second);
        it = matches.erase(it);
      }
      for (const auto& reference : matches)
        orphaned.erase({reference.run_uuid, data.findRunByUuid(reference.run_uuid)->getIdentificationForMatch(reference.match).getId()});
      feature.getIDQueries().insert(orphaned.begin(), orphaned.end());
    }

    /// Remove the links of @p feature to erased identifications and matches (see eraseIdentifications())
    void unlinkErased(BaseFeature& feature, const std::pair<std::set<IdentificationData::QueryReference>, std::set<IdentificationData::MatchReference>>& erased)
    {
      std::erase_if(feature.getIDQueries(), [&](const auto& reference) { return erased.first.contains(reference); });
      std::erase_if(feature.getIDMatches(), [&](const auto& reference) { return erased.second.contains(reference); });
    }
  } // namespace

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
    states(static_cast<size_t>(BaseFeature::AnnotationState::SIZE_OF_ANNOTATIONSTATE), 0) // initialize all with 0
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

  FeatureMap::FeatureMap(const FeatureMap& source):
      MetaInfoInterface(source),
      RangeManagerContainerType(source),
      DocumentIdentifier(source),
      ExposedVector<Feature>(source),
      UniqueIdInterface(source),
      UniqueIdIndexer<FeatureMap>(source),
      protein_identifications_(source.protein_identifications_),
      unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
      data_processing_(source.data_processing_),
      id_data_(source.id_data_)
  {
  }

  FeatureMap::FeatureMap(FeatureMap&& source) = default;

  FeatureMap::~FeatureMap() = default;

  FeatureMap& FeatureMap::operator=(const FeatureMap& rhs) // TODO: cannot be defaulted since OpenMS::IdentificationData is missing operator=
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

    id_data_ = rhs.id_data_;

    return *this;
  }

  FeatureMap& FeatureMap::operator=(FeatureMap&&) = default;


  bool FeatureMap::operator==(const FeatureMap& rhs) const
  {
    return data_ == rhs.data_ && MetaInfoInterface::operator==(rhs) && RangeManagerType::operator==(rhs) && DocumentIdentifier::operator==(rhs)
           && UniqueIdInterface::operator==(rhs) && protein_identifications_ == rhs.protein_identifications_
           && unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ && data_processing_ == rhs.data_processing_
           && id_data_ == rhs.id_data_;
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
    // Check identification compatibility before changing measurements or annotations.
    id_data_.merge(rhs.id_data_);

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

    // append features:
    this->insert(this->end(), rhs.begin(), rhs.end());

    // todo: check for double entries
    // features, unassignedpeptides, proteins...


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
    #ifdef OPENMS_ASSERTIONS
      double rt_min = RangeRT::isEmpty() ? 0 : getMinRT();
      double rt_max = RangeRT::isEmpty() ? 0 : getMaxRT();
      double mz_min = RangeMZ::isEmpty() ? 0 : getMinMZ();
      double mz_max = RangeMZ::isEmpty() ? 0 : getMaxMZ();
      double int_min = RangeIntensity::isEmpty() ? 0 : getMinIntensity();
      double int_max = RangeIntensity::isEmpty() ? 0 : getMaxIntensity();
    #endif

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

    #ifdef OPENMS_ASSERTIONS
      // check if updateRanges() was necessary and find places where it was not
      double rt_min_new = RangeRT::isEmpty() ? 0 : getMinRT();
      double rt_max_new = RangeRT::isEmpty() ? 0 : getMaxRT();
      double mz_min_new = RangeMZ::isEmpty() ? 0 : getMinMZ();
      double mz_max_new = RangeMZ::isEmpty() ? 0 : getMaxMZ();
      double int_min_new = RangeIntensity::isEmpty() ? 0 : getMinIntensity();
      double int_max_new = RangeIntensity::isEmpty() ? 0 : getMaxIntensity();

      // check if all are equal and no update range was necessary
      if (rt_min_new == rt_min && rt_max_new == rt_max
        && int_min_new == int_min && int_max_new == int_max
        && mz_min_new == mz_min && mz_max_new == mz_max)
      {
        OPENMS_LOG_WARN << "Update ranges was called but ranges were already up-to-date" << std::endl;
      }
    #endif
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

  const ProteinIdentification* FeatureMap::findProteinIdentification(const std::string& identifier) const
  {
    for (const auto& prot_id : protein_identifications_)
    {
      if (prot_id.getIdentifier() == identifier)
      {
        return &prot_id;
      }
    }
    return nullptr;
  }

  ProteinIdentification* FeatureMap::findProteinIdentification(const std::string& identifier)
  {
    for (auto& prot_id : protein_identifications_)
    {
      if (prot_id.getIdentifier() == identifier)
      {
        return &prot_id;
      }
    }
    return nullptr;
  }

  const PeptideIdentificationList& FeatureMap::getUnassignedPeptideIdentifications() const
  {
    return unassigned_peptide_identifications_;
  }

  PeptideIdentificationList& FeatureMap::getUnassignedPeptideIdentifications()
  {
    return unassigned_peptide_identifications_;
  }

  void FeatureMap::setUnassignedPeptideIdentifications(const PeptideIdentificationList& unassigned_peptide_identifications)
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

    for (const std::string& filename : s)
    {
      if (!StringUtils::hasSuffix(filename, "mzML") && !StringUtils::hasSuffix(filename, "mzml"))
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
    if (ms_path.size() == 1 && StringUtils::hasSuffix(ms_path[0], "mzML") && File::exists(ms_path[0]))
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
      result += iter->getAnnotationState(id_data_);
    }
    return result;
  }


  std::set<IdentificationData::MatchReference> FeatureMap::getUnassignedIDMatches() const
  {
    std::set<IdentificationData::MatchReference> all, assigned, result;
    for (const auto& run : id_data_.getRuns())
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            all.insert({run.getUuid(), match.getId()});
    const auto collect = [&](const auto& self, const Feature& feature) -> void {
      assigned.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
      for (const auto& subordinate : feature.getSubordinates())
        self(self, subordinate);
    };
    for (const auto& feature : *this)
      collect(collect, feature);
    std::set_difference(all.begin(), all.end(), assigned.begin(), assigned.end(), std::inserter(result, result.end()));
    return result;
  }

  Size FeatureMap::eraseMatches(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&,
                                                        const IdentificationData::Match&)>& remove)
  {
    const auto erased = id_data_.eraseMatches(remove);
    if (erased.empty()) return 0;
    const auto relink = [&](const auto& self, Feature& feature) -> void {
      relinkErased(feature, id_data_, erased);
      for (auto& subordinate : feature.getSubordinates())
        self(self, subordinate);
    };
    for (auto& feature : *this)
      relink(relink, feature);
    return erased.size();
  }

  Size FeatureMap::eraseIdentifications(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&)>& remove)
  {
    const auto erased = id_data_.eraseIdentifications(remove);
    if (erased.first.empty()) return 0;
    const auto unlink = [&](const auto& self, Feature& feature) -> void {
      unlinkErased(feature, erased);
      for (auto& subordinate : feature.getSubordinates())
        self(self, subordinate);
    };
    for (auto& feature : *this)
      unlink(unlink, feature);
    return erased.first.size();
  }


  namespace
  {
    /**
      @brief Erase the features of @p map that @p erase flags, and from its identification data what only they link

      @p links collects the query and match links of a feature (with its subordinates).
    */
    template<class Map, class Links>
    Size eraseFlaggedFeatures(Map& map, const std::vector<bool>& erase, const Links& links)
    {
      const Size n_erased = std::count(erase.begin(), erase.end(), true);
      if (n_erased == 0) return 0;
      auto& data = map.getIdentificationData();
      if (! data.empty())
      {
        std::set<IdentificationData::QueryReference> kept_queries, erased_queries;
        std::set<IdentificationData::MatchReference> kept_matches, erased_matches;
        for (Size i = 0; i < map.size(); ++i)
        {
          if (erase[i]) links(map[i], erased_queries, erased_matches);
          else links(map[i], kept_queries, kept_matches);
        }
        for (const auto& match : kept_matches) erased_matches.erase(match);
        for (const auto& query : kept_queries) erased_queries.erase(query);
        data.eraseMatches([&](const IdentificationData::Run& run, const IdentificationData::Identification&, const IdentificationData::Match& match) {
          return erased_matches.contains({run.getUuid(), match.getId()});
        });
        data.eraseIdentifications([&](const IdentificationData::Run& run, const IdentificationData::Identification& query) {
          return query.getMatches().empty() && erased_queries.contains({run.getUuid(), query.getId()});
        });
      }
      Size index = 0;
      map.erase(std::remove_if(map.begin(), map.end(), [&](const auto&) { return erase[index++]; }), map.end());
      return n_erased;
    }
  } // namespace

  Size FeatureMap::eraseFeatures(const std::function<bool(const Feature&)>& remove)
  {
    std::vector<bool> erase;
    erase.reserve(size());
    for (const auto& feature : *this) erase.push_back(remove(feature));
    const auto links = [&](const Feature& feature, std::set<IdentificationData::QueryReference>& queries,
                           std::set<IdentificationData::MatchReference>& matches) {
      const auto collect = [&](const auto& self, const Feature& item) -> void {
        const auto linked = item.getLinkedIDQueries(id_data_);
        queries.insert(linked.begin(), linked.end());
        matches.insert(item.getIDMatches().begin(), item.getIDMatches().end());
        for (const auto& subordinate : item.getSubordinates()) self(self, subordinate);
      };
      collect(collect, feature);
    };
    return eraseFlaggedFeatures(*this, erase, links);
  }

  Size FeatureMap::eraseUnassignedIdentifications()
  {
    using ID = IdentificationData;
    std::set<ID::QueryReference> linked_queries;
    std::set<ID::MatchReference> linked_matches;
    const auto collect = [&](const auto& self, const Feature& feature) -> void {
      linked_queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
      linked_matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
      for (const auto& subordinate : feature.getSubordinates())
        self(self, subordinate);
    };
    for (const auto& feature : *this)
      collect(collect, feature);
    std::set<ID::QueryReference> queries;
    std::set<ID::MatchReference> matches;
    for (const auto& entry : id_data_.getUnlinked(linked_queries, linked_matches))
    {
      const ID::QueryReference reference {entry.run->getUuid(), entry.query->getId()};
      // An identification that no feature links, with no linked match, is unassigned as a whole. Of one that a feature
      // links (as a peptide identification without hits), or links through some of its matches, the other matches are.
      if (! linked_queries.contains(reference) && entry.matches.size() == entry.query->getMatches().size()) queries.insert(reference);
      else
        for (const auto* match : entry.matches) matches.insert({reference.run_uuid, match->getId()});
    }
    const Size n_matches = eraseMatches([&](const ID::Run& run, const ID::Identification&, const ID::Match& match) {
      return matches.contains({run.getUuid(), match.getId()});
    });
    return n_matches + eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) {
      return queries.contains({run.getUuid(), query.getId()});
    });
  }

  std::vector<IdentificationData::QueryMatches> FeatureMap::getUnassignedIdentifications() const
  {
    std::set<IdentificationData::QueryReference> queries;
    std::set<IdentificationData::MatchReference> matches;
    const auto collect = [&](const auto& self, const Feature& feature) -> void {
      queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
      matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
      for (const auto& subordinate : feature.getSubordinates())
        self(self, subordinate);
    };
    for (const auto& feature : *this)
      collect(collect, feature);
    return id_data_.getUnlinked(queries, matches);
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
