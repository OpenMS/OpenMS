// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>

namespace OpenMS
{
  ConsensusMap::ConsensusMap() = default;

  ConsensusMap::ConsensusMap(const ConsensusMap& source) :
    MetaInfoInterface(source),
    RangeManagerContainerType(source),
    DocumentIdentifier(source),
    ExposedVector<ConsensusFeature>(source),
    UniqueIdInterface(source),
    UniqueIdIndexer<ConsensusMap>(source),
    column_description_(source.column_description_),
    experiment_type_(source.experiment_type_),
    protein_identifications_(source.protein_identifications_),
    unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
    data_processing_(source.data_processing_),
    id_data_() // updated below
  {
    // copy ID data and update references in features:
    IdentificationData::RefTranslator trans = id_data_.merge(source.id_data_);
    for (ConsensusFeature& feature : *this)
    {
      feature.updateIDReferences(trans);
    }
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

    // append proteinIdentification
    protein_identifications_.insert(protein_identifications_.end(),
                                    rhs.protein_identifications_.begin(),
                                    rhs.protein_identifications_.end());

    // ensure non-redundant modification parameter
    for (auto& pi : protein_identifications_)
    {
      std::vector<String>::iterator it_2;

      // remove redundant variable modifications
      std::vector<String>& varMod = pi.getSearchParameters().variable_modifications;
      sort(varMod.begin(), varMod.end());
      it_2 = unique(varMod.begin(), varMod.end());
      varMod.resize(it_2 - varMod.begin());

      // remove redundant fixed modifications
      std::vector<String>& fixMod = pi.getSearchParameters().fixed_modifications;
      sort(fixMod.begin(), fixMod.end());
      it_2 = unique(fixMod.begin(), fixMod.end());
      fixMod.resize(it_2 - fixMod.begin());
    }

    // append unassigned PeptideIdentifications
    unassigned_peptide_identifications_.insert(unassigned_peptide_identifications_.end(),
                                               rhs.unassigned_peptide_identifications_.begin(),
                                               rhs.unassigned_peptide_identifications_.end());

    Size old_size = size();
    // append consensusElements to consensusElementList:
    this->insert(this->end(), rhs.begin(), rhs.end());

    // combine IDs (new format):
    IdentificationData::RefTranslator trans = id_data_.merge(rhs.id_data_);
    // update IDs in new consensus features:
    for (Size i = old_size; i < size(); ++i) {
      (*this)[i].updateIDReferences(trans);
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

  ConsensusMap& ConsensusMap::appendColumns(const ConsensusMap& rhs)
  {
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

    // append proteinIdentification
    protein_identifications_.insert(protein_identifications_.end(),
                                    rhs.protein_identifications_.begin(),
                                    rhs.protein_identifications_.end());

    // ensure non-redundant modification parameter
    for (auto& pi : protein_identifications_)
    {
      std::vector<String>::iterator it_2;

      // remove redundant variable modifications
      std::vector<String>& varMod = pi.getSearchParameters().variable_modifications;
      sort(varMod.begin(), varMod.end());
      it_2 = unique(varMod.begin(), varMod.end());
      varMod.resize(it_2 - varMod.begin());

      // remove redundant fixed modifications
      std::vector<String>& fixMod = pi.getSearchParameters().fixed_modifications;
      sort(fixMod.begin(), fixMod.end());
      it_2 = unique(fixMod.begin(), fixMod.end());
      fixMod.resize(it_2 - fixMod.begin());
    }

    // append unassigned identifications and update map index:
    for (PeptideIdentification pid : rhs.unassigned_peptide_identifications_)
    {
      if (pid.metaValueExists("map_index"))
      {
        Size old_index = pid.getMetaValue("map_index");
        pid.setMetaValue("map_index", lhs_map_size + old_index);
      }
      unassigned_peptide_identifications_.push_back(pid);
    }

    // combine IDs (new format):
    IdentificationData::RefTranslator trans = id_data_.merge(rhs.id_data_);

    // append consensusElements to consensusElementList and update map index:
    for (ConsensusFeature cf : rhs)
    {
      for (PeptideIdentification& pid : cf.getPeptideIdentifications())
      {
        if (pid.metaValueExists("map_index"))
        {
          Size old_index = pid.getMetaValue("map_index");
          pid.setMetaValue("map_index", lhs_map_size + old_index);
        }
      }

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

      // update IDs (new format):
      cf.updateIDReferences(trans);

      emplace_back(cf);
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
      protein_identifications_.clear();
      unassigned_peptide_identifications_.clear();
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

  const String& ConsensusMap::getExperimentType() const
  {
    return experiment_type_;
  }

  void ConsensusMap::setExperimentType(const String& experiment_type)
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

  void ConsensusMap::sortPeptideIdentificationsByMapIndex()
  {
    // lambda predicate
    auto mapIndexLess = [] (const PeptideIdentification & a, const PeptideIdentification & b) -> bool
    {
      const bool has_a = a.metaValueExists("map_index");
      const bool has_b = b.metaValueExists("map_index");

      // moves IDs without meta value to end
      if (has_a && !has_b)
      {
        return true;
      }
      if (!has_a && has_b)
      {
        return false;
      }

      // both have map index annotated
      if (has_a && has_b)
      {
        return a.getMetaValue("map_index") < b.getMetaValue("map_index");
      }

      // no map index annotated in both
      return false;
    };

    std::transform(begin(), end(), begin(),
      [mapIndexLess](ConsensusFeature& c)
      {
        auto& pids = c.getPeptideIdentifications();
        stable_sort(pids.begin(), pids.end(), mapIndexLess);
        return c;
      });
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
    protein_identifications_.swap(from.protein_identifications_);
    unassigned_peptide_identifications_.swap(from.unassigned_peptide_identifications_);
    data_processing_.swap(from.data_processing_);
    id_data_.swap(from.id_data_);
  }

  /// non-mutable access to the protein identifications
  const std::vector<ProteinIdentification>& ConsensusMap::getProteinIdentifications() const
  {
    return protein_identifications_;
  }

  /// mutable access to the protein identifications
  std::vector<ProteinIdentification>& ConsensusMap::getProteinIdentifications()
  {
    return protein_identifications_;
  }

  /// sets the protein identifications
  void ConsensusMap::setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)
  {
    protein_identifications_ = protein_identifications;
  }

  /// sets the protein identifications
  void ConsensusMap::setProteinIdentifications(std::vector<ProteinIdentification>&& protein_identifications)
  {
    protein_identifications_ = std::move(protein_identifications);
  }

  /// non-mutable access to the unassigned peptide identifications
  const std::vector<PeptideIdentification>& ConsensusMap::getUnassignedPeptideIdentifications() const
  {
    return unassigned_peptide_identifications_;
  }

  /// mutable access to the unassigned peptide identifications
  std::vector<PeptideIdentification>& ConsensusMap::getUnassignedPeptideIdentifications()
  {
    return unassigned_peptide_identifications_;
  }

  /// sets the unassigned peptide identifications
  void ConsensusMap::setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications)
  {
    unassigned_peptide_identifications_ = unassigned_peptide_identifications;
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
      OPENMS_LOG_WARN << "Setting empty MS runs paths. Expected one for each map. Resulting ConsensusMap contains " + String(column_description_.size()) + " maps." << std::endl;
      for (auto & cd : column_description_)
      {
        cd.second.filename = "UNKNOWN";
      }
    }
    else if (!column_description_.empty() && s.size() != column_description_.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Number of MS runs paths (" + String(s.size()) +
        ") must match number of columns (" + String(column_description_.size()) + ").");
    }

    Size i(0);
    for (auto const & p : s)
    {
      if (!p.hasSuffix("mzML") && !p.hasSuffix("mzml"))
      {
        OPENMS_LOG_WARN << "To ensure tracability of results please prefer mzML files as primary MS run." << std::endl
                        << "Filename: '" << p << "'" << std::endl;
      }

      column_description_[i].filename = p;
      ++i;
    }
  }

  void ConsensusMap::setPrimaryMSRunPath(const StringList& s, MSExperiment & e)
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
    return data_ == rhs.data_ &&
           MetaInfoInterface::operator==(rhs) &&
           RangeManagerType::operator==(rhs) &&
           DocumentIdentifier::operator==(rhs) &&
           UniqueIdInterface::operator==(rhs) &&
           column_description_ == rhs.column_description_ &&
           experiment_type_ == rhs.experiment_type_ &&
           protein_identifications_ == rhs.protein_identifications_ &&
           unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ &&
           data_processing_ == rhs.data_processing_;
    // @TODO: implement "operator==" for IdentificationData?
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
  }

  bool ConsensusMap::isMapConsistent(Logger::LogStream* stream) const
  {
    Size stats_wrongMID(0); // invalid map ID references by a feature handle
    std::map<Size, Size> wrong_ID_count; // which IDs were given which are not valid

    // check file descriptions
    std::set<String> maps;
    String all_maps; // for output later
    for (ColumnHeaders::const_iterator it = column_description_.begin(); it != column_description_.end(); ++it)
    {
      String s = String("  file: ") + it->second.filename + " label: " + it->second.label;
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
        if (column_description_.find(it->getMapIndex()) == column_description_.end())
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
    // @TODO: handle IDs in new format (IdentificationData)

    Size numbr_exps = column_description_.size();
    std::vector<FeatureMap>fmaps(numbr_exps);

    // Check for Isobaric Analyzer
    bool iso_analyze = QCBase::isLabeledExperiment(*this);

    for (const auto& cf : *this)
    {
      UInt64 min_index = std::numeric_limits<UInt64>::max();
      // Create new Features from FeatureHandles
      std::map<UInt64, BaseFeature> new_feats;
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
      for (const PeptideIdentification& pep_id : cf.getPeptideIdentifications())
      {
        // ... the first Feature.
        if (iso_analyze)
        {
          (*new_feats.begin()).second.getPeptideIdentifications().push_back(pep_id);
          continue;
        }

        // ... the corresponding Feature by map_index.
        if (!pep_id.metaValueExists("map_index"))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "File did not undergo IsobaricAnalyzer, but no map index was found at PeptideIdentifications. Check Input!");
        }
        new_feats[pep_id.getMetaValue("map_index")].getPeptideIdentifications().push_back(pep_id);
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

      // Add new Features to corresponding FeatureMap.
      for (auto it = new_feats.begin(); it != new_feats.end(); ++it)
      {
        fmaps[it->first].emplace_back(std::move(it->second));
      }
    }

    // Add unassigned PeptideIdentifications to ...
    if (iso_analyze)
    {
      // ... the first FeatureMap.
      fmaps[0].getUnassignedPeptideIdentifications() = this->getUnassignedPeptideIdentifications();
      fmaps[0].getProteinIdentifications() = this->getProteinIdentifications(); // wrong! improve: only copy the ProtID which belongs to this FMap!
    }
    else
    {
      // ... the corresponding FeatureMap by map_index.
      for (const PeptideIdentification& upep_id : this->getUnassignedPeptideIdentifications())
      {
        if (!upep_id.metaValueExists("map_index"))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "File did not undergo IsobaricAnalyzer, but no map index was found at PeptideIdentifications. Check Input!");
        }
        fmaps[upep_id.getMetaValue("map_index")].getUnassignedPeptideIdentifications().push_back(upep_id);
      }
    }

    for (auto& fm : fmaps)
    {
      fm.getDataProcessing() = this->getDataProcessing();
    }

    return fmaps;
  }

  unsigned ConsensusMap::ColumnHeader::getLabelAsUInt(const String& experiment_type) const
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

  std::set<IdentificationDataInternal::ObservationMatchRef> ConsensusMap::getUnassignedIDMatches() const
  {
    std::set<IdentificationData::ObservationMatchRef> all_matches;
    for (auto it = id_data_.getObservationMatches().begin();
         it != id_data_.getObservationMatches().end(); ++it)
    {
      all_matches.insert(it);
    }
    std::set<IdentificationData::ObservationMatchRef> assigned_matches;
    for (const ConsensusFeature& feat : *this)
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


  const IdentificationData& ConsensusMap::getIdentificationData() const
  {
    return id_data_;
  }


  IdentificationData& ConsensusMap::getIdentificationData()
  {
    return id_data_;
  }

} // namespace OpenMS
