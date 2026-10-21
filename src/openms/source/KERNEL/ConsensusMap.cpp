// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/DataProcessingUtils.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>

namespace OpenMS
{

  namespace
  {
    /// Sort the fixed and variable modifications of every run and remove duplicates (as merging maps does for their
    /// protein identifications)
    void uniqueModifications(IdentificationData& data)
    {
      std::vector<std::string> identifiers;
      for (const auto& run : data.getRuns()) identifiers.push_back(run.getIdentifier());
      for (const auto& identifier : identifiers)
      {
        auto& run = data.getRun(identifier);
        auto settings = run.getSettings();
        for (auto* mods : {&settings.search.fixed_modifications, &settings.search.variable_modifications})
        {
          std::sort(mods->begin(), mods->end());
          mods->erase(std::unique(mods->begin(), mods->end()), mods->end());
        }
        if (! (settings == run.getSettings())) run.setSettings(settings);
      }
      // (the protein runs of inference results carry search parameters as well, which export writes)
      auto results = data.getInferenceResults();
      bool changed = false;
      for (auto& result : results)
      {
        auto& search = result.proteins.getSearchParameters();
        for (auto* mods : {&search.fixed_modifications, &search.variable_modifications})
        {
          auto unique_mods = *mods;
          std::sort(unique_mods.begin(), unique_mods.end());
          unique_mods.erase(std::unique(unique_mods.begin(), unique_mods.end()), unique_mods.end());
          if (unique_mods != *mods)
          {
            *mods = std::move(unique_mods);
            changed = true;
          }
        }
      }
      if (! changed) return;
      data.clearInferenceResults();
      for (auto& result : results) data.addInferenceResult(std::move(result));
    }

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

  ConsensusMap::ConsensusMap() = default;

  ConsensusMap::ConsensusMap(const ConsensusMap& source):
      MetaInfoInterface(source),
      RangeManagerContainerType(source),
      DocumentIdentifier(source),
      ExposedVector<ConsensusFeature>(source),
      UniqueIdInterface(source),
      UniqueIdIndexer<ConsensusMap>(source),
      column_description_(source.column_description_),
      experiment_type_(source.experiment_type_),
      data_processing_(source.data_processing_),
      id_data_(source.id_data_)
  {
  }

  ConsensusMap::ConsensusMap(ConsensusMap&& source) = default;

  ConsensusMap::~ConsensusMap() = default;

  ConsensusMap::ConsensusMap(size_type n) :
    ExposedVector<ConsensusFeature>(n)
  {
  }

  ConsensusMap& ConsensusMap::operator=(const ConsensusMap& source) = default;

  ConsensusMap& ConsensusMap::appendRows(const ConsensusMap& rhs)
  {
    // Check identification compatibility before changing measurements or annotations.
    id_data_.merge(rhs.id_data_);

    ConsensusMap empty_map;

    // reset these:
    RangeManagerContainerType::operator=(empty_map);

    if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty())
    {
      OPENMS_LOG_INFO << "DocumentIdentifiers are lost during merge of ConsensusMaps\n";
    }

    DocumentIdentifier::operator=(empty_map);
    UniqueIdInterface::operator=(empty_map);

    // append dataProcessing
    data_processing_.insert(data_processing_.end(),
                            rhs.data_processing_.begin(),
                            rhs.data_processing_.end());

    // append fileDescription
    column_description_.insert(rhs.column_description_.begin(), rhs.column_description_.end());

    // update filename and map size
    std::map<UInt64, ColumnHeader>::const_iterator it = column_description_.begin();
    std::map<UInt64, ColumnHeader>::const_iterator it2 = rhs.column_description_.begin();

    for (; it != column_description_.end() && it2 != rhs.column_description_.end(); ++it, ++it2)
    {
      getColumnHeaders()[it->first].filename = "mergedConsensusXMLFile";
      getColumnHeaders()[it->first].size = it->second.size + it2->second.size;
    }

    // ensure non-redundant modification parameters in the runs of the identification data
    uniqueModifications(id_data_);

    // append consensusElements to consensusElementList:
    this->insert(this->end(), rhs.begin(), rhs.end());


    // consistency
    try
    {
      UniqueIdIndexer<ConsensusMap>::updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition&) // assign new UID's for conflicting entries
    {
      Size replaced_uids =  UniqueIdIndexer<ConsensusMap>::resolveUniqueIdConflicts();
      OPENMS_LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
    }

    return *this;
  }

  ConsensusMap& ConsensusMap::appendColumns(const ConsensusMap& rhs)
  {
    // The identifications of rhs refer to its columns: shift their map index like the column indices below.
    // Merging checks identification compatibility before measurements or annotations change.
    IdentificationData rhs_data = rhs.id_data_;
    for (const auto& current : rhs.id_data_.getRuns())
    {
      auto& run = rhs_data.getRun(current.getIdentifier());
      for (const auto& source : current.getSources())
        for (const auto& query : source.identifications)
        {
          if (! query.metaValueExists("map_index")) continue;
          if (id_data_.findRunByUuid(run.getUuid()))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Both consensus maps contain identification run '" + run.getIdentifier()
                                            + "' with map indices; its identifications cannot refer to columns of both maps.",
                                          run.getUuid());
          }
          IdentificationData::Observation observation = query;
          observation.setMetaValue("map_index", column_description_.size() + static_cast<Size>(query.getMetaValue("map_index")));
          run.replaceObservation(query.getId(), observation);
        }
    }
    id_data_.merge(rhs_data);

    ConsensusMap empty_map;

    // reset these:
    RangeManagerType::operator=(empty_map);

    if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty())
    {
      OPENMS_LOG_INFO << "DocumentIdentifiers are lost during merge of ConsensusMaps\n";
    }

    DocumentIdentifier::operator=(empty_map);
    UniqueIdInterface::operator=(empty_map);

    // append dataProcessing
    data_processing_.insert(data_processing_.end(),
                            rhs.data_processing_.begin(),
                            rhs.data_processing_.end());

    // append column headers (file descriptions) and increase column index (map index)
    Size lhs_map_size = column_description_.size();
    for (const auto& rhsfd : rhs.column_description_)
    {
      column_description_.insert(std::make_pair(lhs_map_size + rhsfd.first, rhsfd.second));
    }

    // ensure non-redundant modification parameters in the runs of the identification data
    uniqueModifications(id_data_);

    // append consensusElements to consensusElementList and update map index:
    for (ConsensusFeature cf : rhs)
    {
      // update map indices
      ConsensusFeature::HandleSetType new_handles;
      for (auto handle : cf) // OMS_CODING_TEST_EXCLUDE Note: std::set only provides const iterators, so we copy
      {
        //since we only add a constant to the map_index, the set order will not change.
        handle.setMapIndex(lhs_map_size + handle.getMapIndex());
        new_handles.insert(handle);
      }
      cf.setFeatures(std::move(new_handles));
      new_handles.clear();

      // match-between-runs transfers (PipEchoAlgorithm) name their target runs by map index, too
      if (cf.metaValueExists(Constants::UserParam::MBR_TRANSFER_MAP_INDICES))
      {
        IntList transfer_map_indices = cf.getMetaValue(Constants::UserParam::MBR_TRANSFER_MAP_INDICES);
        for (Int& map_index : transfer_map_indices)
        {
          map_index += static_cast<Int>(lhs_map_size);
        }
        cf.setMetaValue(Constants::UserParam::MBR_TRANSFER_MAP_INDICES, transfer_map_indices);
      }

      emplace_back(std::move(cf));
    }

    // consistency
    try
    {
      UniqueIdIndexer<ConsensusMap>::updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition&) // assign new UID's for conflicting entries
    {
      Size replaced_uids =  UniqueIdIndexer<ConsensusMap>::resolveUniqueIdConflicts();
      OPENMS_LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
    }

    return *this;
  }


  void ConsensusMap::clear(bool clear_meta_data)
  {
    data_.clear();

    if (clear_meta_data)
    {
      clearMetaInfo();
      clearRanges();
      // no "clear" method for DocumentIdentifier available
      this->DocumentIdentifier::operator=(DocumentIdentifier());
      clearUniqueId();
      column_description_.clear();
      experiment_type_ = "label-free";  // default
      data_processing_.clear();
      id_data_.clear();
    }
  }

  const ConsensusMap::ColumnHeaders& ConsensusMap::getColumnHeaders() const
  {
    return column_description_;
  }

  ConsensusMap::ColumnHeaders& ConsensusMap::getColumnHeaders()
  {
    return column_description_;
  }

  void ConsensusMap::setColumnHeaders(const ConsensusMap::ColumnHeaders& column_description)
  {
    column_description_ = column_description;
  }

  const std::string& ConsensusMap::getExperimentType() const
  {
    return experiment_type_;
  }

  void ConsensusMap::setExperimentType(const std::string& experiment_type)
  {
    if (experiment_type != "label-free"
      && experiment_type != "labeled_MS1"
      && experiment_type != "labeled_MS2")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unknown experiment type. " + experiment_type + ". Must be one of (label-free, labeled_MS1, labeled_MS2)");
    }
    experiment_type_ = experiment_type;
  }

  void ConsensusMap::sortByIntensity(bool reverse)
  {
    if (reverse)
    {
      std::stable_sort(begin(), end(), [](auto &left, auto &right) {ConsensusFeature::IntensityLess cmp; return cmp(right, left);});
    }
    else
    {
      std::stable_sort(begin(), end(), ConsensusFeature::IntensityLess());
    }
  }

  void ConsensusMap::sortByRT()
  {
    std::stable_sort(begin(), end(), ConsensusFeature::RTLess());
  }

  void ConsensusMap::sortByMZ()
  {
    std::stable_sort(begin(), end(), ConsensusFeature::MZLess());
  }

  void ConsensusMap::sortByPosition()
  {
    std::stable_sort(begin(), end(), ConsensusFeature::PositionLess());
  }

  void ConsensusMap::sortByQuality(bool reverse)
  {
    if (reverse)
    {
      std::stable_sort(begin(), end(), [](auto &left, auto &right) {ConsensusFeature::QualityLess cmp; return cmp(right, left);});
    }
    else
    {
      std::stable_sort(begin(), end(), ConsensusFeature::QualityLess());
    }
  }

  void ConsensusMap::sortBySize()
  {
    std::stable_sort(begin(), end(), [](auto &left, auto &right) {ConsensusFeature::SizeLess cmp; return cmp(right, left);});
  }

  void ConsensusMap::sortByMaps()
  {
    std::stable_sort(begin(), end(), ConsensusFeature::MapsLess());
  }

  void ConsensusMap::swap(ConsensusMap& from)
  {
    ConsensusMap tmp;

    //swap range information
    tmp.RangeManagerType::operator=(* this);
    this->RangeManagerType::operator=(from);
    from.RangeManagerType::operator=(tmp);

    //swap consensus features
    data_.swap(from.data_);

    // swap DocumentIdentifier
    DocumentIdentifier::swap(from);

    // swap unique id
    UniqueIdInterface::swap(from);

    // swap unique id index
    UniqueIdIndexer<ConsensusMap>::swap(from);

    // swap the remaining members
    std::swap(column_description_, from.column_description_);
    experiment_type_.swap(from.experiment_type_);
    data_processing_.swap(from.data_processing_);
    id_data_.swap(from.id_data_);
  }


  /// returns a const reference to the description of the applied data processing
  const std::vector<DataProcessing>& ConsensusMap::getDataProcessing() const
  {
    return data_processing_;
  }

  /// returns a mutable reference to the description of the applied data processing
  std::vector<DataProcessing>& ConsensusMap::getDataProcessing()
  {
    return data_processing_;
  }

  /// sets the description of the applied data processing
  void ConsensusMap::setDataProcessing(const std::vector<DataProcessing>& processing_method)
  {
    data_processing_ = processing_method;
  }

  /// set the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
  void ConsensusMap::setPrimaryMSRunPath(const StringList& s)
  {
    if (s.empty())
    {
      OPENMS_LOG_WARN << "Setting empty MS runs paths. Expected one for each map. Resulting ConsensusMap contains " + StringUtils::toStr(column_description_.size()) + " maps." << std::endl;
      for (auto & cd : column_description_)
      {
        cd.second.filename = "UNKNOWN";
      }
    }
    else if (!column_description_.empty() && s.size() != column_description_.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Number of MS runs paths (" + StringUtils::toStr(s.size()) +
        ") must match number of columns (" + StringUtils::toStr(column_description_.size()) + ").");
    }

    // column_description_ is keyed by map index, and the keys need not be contiguous (e.g.
    // after 'FileFilter -consensus:map 0 3'). When columns already exist, column_description_[i]
    // would silently create a new entry at key i instead of updating the real i-th column by key
    // order, growing the map and leaving the intended column's filename unset.
    bool had_existing_columns = !column_description_.empty();
    auto it = column_description_.begin();
    Size i(0);
    for (auto const & p : s)
    {
      if (!StringUtils::hasSuffix(p, "mzML") && !StringUtils::hasSuffix(p, "mzml"))
      {
        OPENMS_LOG_WARN << "To ensure tracability of results please prefer mzML files as primary MS run." << std::endl
                        << "Filename: '" << p << "'" << std::endl;
      }

      if (had_existing_columns)
      {
        it->second.filename = p;
        ++it;
      }
      else
      {
        column_description_[i].filename = p;
      }
      ++i;
    }
  }

  void ConsensusMap::setPrimaryMSRunPath(const StringList& s, MSExperiment & e)
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

  void ConsensusMap::getPrimaryMSRunPath(StringList& toFill) const
  {
    /// get the file path to the MS run
    for (auto const & fd : column_description_)
    {
      toFill.push_back(fd.second.filename);
    }
  }

  /// Equality operator
  bool ConsensusMap::operator==(const ConsensusMap& rhs) const
  {
    return data_ == rhs.data_ && MetaInfoInterface::operator==(rhs) && RangeManagerType::operator==(rhs) && DocumentIdentifier::operator==(rhs)
           && UniqueIdInterface::operator==(rhs) && column_description_ == rhs.column_description_ && experiment_type_ == rhs.experiment_type_
           && data_processing_ == rhs.data_processing_ && id_data_ == rhs.id_data_;
  }

  /// Equality operator
  bool ConsensusMap::operator!=(const ConsensusMap& rhs) const
  {
    return !(operator==(rhs));
  }

  std::ostream& operator<<(std::ostream& os, const ConsensusMap& cons_map)
  {
    for (ConsensusMap::ColumnHeaders::const_iterator it = cons_map.getColumnHeaders().begin(); it != cons_map.getColumnHeaders().end(); ++it)
    {
      os << "Map " << it->first << ": " << it->second.filename << " - " << it->second.label << " - " << it->second.size << std::endl;
    }

    for (Size i = 0; i < cons_map.size(); ++i)
    {
      os << cons_map[i] << std::endl;
    }

    return os;
  }

  void ConsensusMap::updateRanges()
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
    // enlarge the range by the internal points of each feature
    for (const auto& cf : data_)
    {
      extendRT(cf.getRT());
      extendMZ(cf.getMZ());
      extendIntensity(cf.getIntensity());
      for (const auto& handle : cf.getFeatures())
      {
        extendRT(handle.getRT());
        extendMZ(handle.getMZ());
        extendIntensity(handle.getIntensity());
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

  bool ConsensusMap::isMapConsistent(Logger::LogStream* stream) const
  {
    Size stats_wrongMID(0); // invalid map ID references by a feature handle
    std::map<Size, Size> wrong_ID_count; // which IDs were given which are not valid

    // check file descriptions
    std::set<std::string> maps;
    std::string all_maps; // for output later
    for (ColumnHeaders::const_iterator it = column_description_.begin(); it != column_description_.end(); ++it)
    {
      std::string s =std::string("  file: ") + it->second.filename + " label: " + it->second.label;
      maps.insert(s);
      all_maps += s;
    }

    if (maps.size() != column_description_.size())
    {
      if (stream != nullptr)
      {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
        *stream << "Map descriptions (file name + label) in ConsensusMap are not unique:\n" << all_maps << std::endl;
      }
      return false;
    }

    // check map IDs
    for (Size i = 0; i < size(); ++i)
    {
      const ConsensusFeature& elem = (*this)[i];
      for (ConsensusFeature::HandleSetType::const_iterator it = elem.begin(); it != elem.end(); ++it)
      {
        if (!column_description_.contains(it->getMapIndex()))
        {
          ++stats_wrongMID;
          ++wrong_ID_count[it->getMapIndex()];
        }
      }
    }

    if (stats_wrongMID > 0)
    {
      if (stream != nullptr)
      {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
        *stream << "ConsensusMap contains " << stats_wrongMID << " invalid references to maps:\n";
        for (std::map<Size, Size>::const_iterator it = wrong_ID_count.begin(); it != wrong_ID_count.end(); ++it)
        {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
          *stream << "  wrong id=" << it->first << " (occurred " << it->second << "x)\n";
        }
OPENMS_THREAD_CRITICAL(LOGSTREAM)
        *stream << std::endl;
      }
      return false;
    }

    return true;
  }

  std::vector<FeatureMap> ConsensusMap::split(ConsensusMap::SplitMeta mode) const
  {
    // Column headers are keyed by map index, and the keys need not be contiguous (e.g. after
    // 'FileFilter -consensus:map 0 3', or for a consensusXML without <map id="0">). The result
    // vector is therefore addressed by the position of a map index among the column headers
    // (in key order), never by the raw map index, which may lie outside the vector.
    std::map<UInt64, Size> index_to_position;
    Size position = 0;
    for (const auto& column : column_description_)
    {
      index_to_position[column.first] = position++;
    }
    auto positionOf = [&index_to_position](UInt64 map_index) -> Size
    {
      const auto pos = index_to_position.find(map_index);
      if (pos == index_to_position.end())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Map index " + std::to_string(map_index) + " does not name a column of this ConsensusMap. Check input!");
      }
      return pos->second;
    };
    auto mapIndexOf = [](const PeptideIdentification& pep_id) -> UInt64
    {
      if (!pep_id.metaValueExists("map_index"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "File did not undergo IsobaricAnalyzer, but no map index was found at PeptideIdentifications. Check Input!");
      }
      return static_cast<UInt64>(pep_id.getMetaValue("map_index"));
    };

    std::vector<FeatureMap> fmaps(index_to_position.size());
    // The identifications are split as peptide identifications; each FeatureMap imports its share at the end.
    const MapPeptideIdentifications identifications = IdentificationDataConverter::exportMapIdentifications(*this);
    std::vector<MapPeptideIdentifications> parts(fmaps.size());

    // Check for Isobaric Analyzer
    bool iso_analyze = DataProcessingUtils::hasIsobaricAnalyzer(getDataProcessing());

    for (Size cf_index = 0; cf_index < size(); ++cf_index)
    {
      const ConsensusFeature& cf = (*this)[cf_index];
      UInt64 min_index = std::numeric_limits<UInt64>::max();
      // Create new Features from FeatureHandles
      std::map<UInt64, BaseFeature> new_feats;
      std::map<UInt64, PeptideIdentificationList> new_peptides;
      for (const FeatureHandle& fh : cf.getFeatures())
      {
        UInt64 index = fh.getMapIndex();
        // GCC-OPT 4.8 does not compile with:  new_feats.emplace(index, fh);
        // , thus we use:
        new_feats[index] = BaseFeature(fh);
        min_index = std::min(index, min_index);
      }

      if (iso_analyze)
      {
        if (min_index != 0)
        {
          throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "File seems to have gone through IsobaricAnalyzer, but there was no feature with map index 0 found. Check Input!");
        }
      }

      // Add PeptideIdentifications to ...
      for (const PeptideIdentification& pep_id : identifications.at({cf_index}))
      {
        // ... the first Feature, or the corresponding Feature by map_index.
        const UInt64 index = iso_analyze ? new_feats.begin()->first : mapIndexOf(pep_id);
        if (!new_feats.contains(index)) new_feats[index] = BaseFeature();
        new_peptides[index].push_back(pep_id);
      }

      // handle MetaValues of current CF
      switch (mode)
      {
        case SplitMeta::DISCARD :
          break;

        case SplitMeta::COPY_ALL :
          for (auto it = new_feats.begin(); it != new_feats.end(); ++it)
          {
            (it->second).MetaInfoInterface::operator=(cf);
          }
          break;

        case SplitMeta::COPY_FIRST :
          if (min_index != 0)
          {
            throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "No feature with map index 0 to copy MetaValues to. Check Input or switch mode!");
          }
          new_feats.begin()->second.MetaInfoInterface::operator=(cf);
          break;
      }

      // Add new Features (and their identifications) to corresponding FeatureMap.
      for (auto it = new_feats.begin(); it != new_feats.end(); ++it)
      {
        const Size part = positionOf(it->first);
        if (auto peptides = new_peptides.find(it->first); peptides != new_peptides.end())
        {
          parts[part].features[{fmaps[part].size()}] = std::move(peptides->second);
        }
        fmaps[part].emplace_back(std::move(it->second));
      }
    }

    // Add unassigned PeptideIdentifications to ...
    if (iso_analyze)
    {
      // ... the first FeatureMap, i.e. the one of map index 0. It is resolved like any other map
      // index, so a missing column 0 (or no column at all) is reported instead of written out of bounds.
      const Size first = positionOf(0);
      parts[first].unassigned = identifications.unassigned;
      parts[first].proteins = identifications.proteins; // wrong! improve: only copy the ProtID which belongs to this FMap!
    }
    else
    {
      // ... the corresponding FeatureMap by map_index.
      for (const PeptideIdentification& upep_id : identifications.unassigned)
      {
        parts[positionOf(mapIndexOf(upep_id))].unassigned.push_back(upep_id);
      }

      // every peptide identification needs its protein identification run: copy the runs each FeatureMap refers to
      for (auto& part : parts)
      {
        std::set<std::string> identifiers;
        for (const PeptideIdentification* pep_id : part.peptides()) identifiers.insert(pep_id->getIdentifier());
        for (const ProteinIdentification& run : identifications.proteins)
        {
          if (identifiers.contains(run.getIdentifier())) part.proteins.push_back(run);
        }
      }
    }

    for (Size i = 0; i < fmaps.size(); ++i)
    {
      if (!parts[i].empty()) IdentificationDataConverter::importMapIdentifications(fmaps[i], parts[i]);
      fmaps[i].getDataProcessing() = this->getDataProcessing();
    }

    return fmaps;
  }

  unsigned ConsensusMap::ColumnHeader::getLabelAsUInt(const std::string& experiment_type) const
  {
    if (metaValueExists("channel_id"))
    {
      return static_cast<unsigned int>(getMetaValue("channel_id")) + 1;
    }
    else
    {
      if (experiment_type != "label-free")
      {
        // TODO There seem to be files in our test data from the Multiplex toolset that do not annotate
        //  a channel id but only add the "label" attribute with the SILAC modification. Add a fall-back here?
        OPENMS_LOG_WARN << "No channel id annotated in labelled consensusXML. Assuming only a single channel was used." << std::endl;
      }
      return 1;
    }
  }

  std::set<IdentificationData::MatchReference> ConsensusMap::getUnassignedIDMatches() const
  {
    std::set<IdentificationData::MatchReference> all, assigned, result;
    for (const auto& run : id_data_.getRuns())
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            all.insert({run.getUuid(), match.getId()});
    for (const auto& feature : *this)
      assigned.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
    std::set_difference(all.begin(), all.end(), assigned.begin(), assigned.end(), std::inserter(result, result.end()));
    return result;
  }

  Size ConsensusMap::eraseMatches(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&,
                                                          const IdentificationData::Match&)>& remove)
  {
    const auto erased = id_data_.eraseMatches(remove);
    for (auto& feature : *this)
      if (! erased.empty()) relinkErased(feature, id_data_, erased);
    return erased.size();
  }

  Size ConsensusMap::eraseIdentifications(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&)>& remove)
  {
    const auto erased = id_data_.eraseIdentifications(remove);
    if (erased.first.empty()) return 0;
    for (auto& feature : *this)
      unlinkErased(feature, erased);
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

  Size ConsensusMap::eraseFeatures(const std::function<bool(const ConsensusFeature&)>& remove)
  {
    std::vector<bool> erase;
    erase.reserve(size());
    for (const auto& feature : *this) erase.push_back(remove(feature));
    const auto links = [&](const ConsensusFeature& feature, std::set<IdentificationData::QueryReference>& queries,
                           std::set<IdentificationData::MatchReference>& matches) {
      const auto linked = feature.getLinkedIDQueries(id_data_);
      queries.insert(linked.begin(), linked.end());
      matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
    };
    return eraseFlaggedFeatures(*this, erase, links);
  }

  void ConsensusMap::clearIdentifications()
  {
    for (auto& feature : *this)
    {
      feature.getIDQueries().clear();
      feature.getIDMatches().clear();
    }
    id_data_.clear();
  }

  Size ConsensusMap::eraseUnassignedIdentifications()
  {
    using ID = IdentificationData;
    std::set<ID::QueryReference> linked_queries;
    std::set<ID::MatchReference> linked_matches;
    for (const auto& feature : *this)
    {
      linked_queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
      linked_matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
    }
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

  std::vector<IdentificationData::QueryMatches> ConsensusMap::getUnassignedIdentifications() const
  {
    std::set<IdentificationData::QueryReference> queries;
    std::set<IdentificationData::MatchReference> matches;
    for (const auto& feature : *this)
    {
      queries.insert(feature.getIDQueries().begin(), feature.getIDQueries().end());
      matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
    }
    return id_data_.getUnlinked(queries, matches);
  }

  const IdentificationData& ConsensusMap::getIdentificationData() const
  {
    return id_data_;
  }


  IdentificationData& ConsensusMap::getIdentificationData()
  {
    return id_data_;
  }

} // namespace OpenMS
