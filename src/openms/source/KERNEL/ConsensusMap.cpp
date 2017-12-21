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

#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{

  ConsensusMap::FileDescription::FileDescription() :
    MetaInfoInterface(),
    filename(),
    label(),
    size(0),
    unique_id(UniqueIdInterface::INVALID)
  {
  }

  ConsensusMap::FileDescription::FileDescription(const ConsensusMap::FileDescription& other) :
    MetaInfoInterface(other),
    filename(other.filename),
    label(other.label),
    size(other.size),
    unique_id(other.unique_id)
  {
  }

  ConsensusMap::ConsensusMap() :
    Base(),
    MetaInfoInterface(),
    RangeManagerType(),
    DocumentIdentifier(),
    UniqueIdInterface(),
    UniqueIdIndexer<ConsensusMap>(),
    file_description_(),
    experiment_type_(),
    protein_identifications_(),
    unassigned_peptide_identifications_(),
    data_processing_()
  {
  }

  ConsensusMap::ConsensusMap(const ConsensusMap& source) :
    Base(source),
    MetaInfoInterface(source),
    RangeManagerType(source),
    DocumentIdentifier(source),
    UniqueIdInterface(source),
    UniqueIdIndexer<ConsensusMap>(source),
    file_description_(source.file_description_),
    experiment_type_(source.experiment_type_),
    protein_identifications_(source.protein_identifications_),
    unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
    data_processing_(source.data_processing_)
  {
  }

  ConsensusMap::~ConsensusMap()
  {
  }

  ConsensusMap::ConsensusMap(Base::size_type n) :
    Base(n),
    MetaInfoInterface(),
    RangeManagerType(),
    DocumentIdentifier(),
    UniqueIdInterface(),
    file_description_(),
    experiment_type_(),
    protein_identifications_(),
    unassigned_peptide_identifications_(),
    data_processing_()
  {
  }

  ConsensusMap& ConsensusMap::operator=(const ConsensusMap& source)
  {
    if (this == &source)
    {
      return *this;
    }

    Base::operator=(source);
    MetaInfoInterface::operator=(source);
    RangeManagerType::operator=(source);
    DocumentIdentifier::operator=(source);
    UniqueIdInterface::operator=(source);
    file_description_ = source.file_description_;
    experiment_type_ = source.experiment_type_;
    protein_identifications_ = source.protein_identifications_;
    unassigned_peptide_identifications_ = source.unassigned_peptide_identifications_;
    data_processing_ = source.data_processing_;

    return *this;
  }

  ConsensusMap& ConsensusMap::operator+=(const ConsensusMap& rhs)
  {
    ConsensusMap empty_map;

    // reset these:
    RangeManagerType::operator=(empty_map);

    if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty())
    {
      LOG_INFO << "DocumentIdentifiers are lost during merge of ConsensusMaps\n";
    }

    DocumentIdentifier::operator=(empty_map);
    UniqueIdInterface::operator=(empty_map);

    // append spectra_data information
    StringList thisRuns_;
    this->getPrimaryMSRunPath(thisRuns_);
    StringList rhsRuns_;
    rhs.getPrimaryMSRunPath(rhsRuns_);
    thisRuns_.insert(thisRuns_.end(), rhsRuns_.begin(), rhsRuns_.end());
    this->setPrimaryMSRunPath(thisRuns_);

    // append dataProcessing
    data_processing_.insert(data_processing_.end(),
                            rhs.data_processing_.begin(),
                            rhs.data_processing_.end());

    // append fileDescription
    file_description_.insert(rhs.file_description_.begin(), rhs.file_description_.end());

    // update filename and map size
    Map<UInt64, FileDescription>::const_iterator it = file_description_.begin();
    Map<UInt64, FileDescription>::const_iterator it2 = rhs.file_description_.begin();

    for (; it != file_description_.end() && it2 != rhs.file_description_.end(); ++it, ++it2)
    {
      getFileDescriptions()[it->first].filename = "mergedConsensusXMLFile";
      getFileDescriptions()[it->first].size = it->second.size + it2->second.size;
    }

    // append proteinIdentification
    protein_identifications_.insert(protein_identifications_.end(),
                                    rhs.protein_identifications_.begin(),
                                    rhs.protein_identifications_.end());

    // ensure non-redundant modification parameter
    for (std::vector<ProteinIdentification>::iterator it_1 = protein_identifications_.begin();
         it_1 != protein_identifications_.end();
         ++it_1)
    {
      std::vector<String>::iterator it_2;

      // remove redundant variable modifications
      std::vector<String>& varMod = const_cast<std::vector<String>&>(it_1->getSearchParameters().variable_modifications);
      sort(varMod.begin(), varMod.end());
      it_2 = unique(varMod.begin(), varMod.end());
      varMod.resize(it_2 - varMod.begin());

      // remove redundant fixed modifications
      std::vector<String>& fixMod = const_cast<std::vector<String>&>(it_1->getSearchParameters().fixed_modifications);
      sort(fixMod.begin(), fixMod.end());
      it_2 = unique(fixMod.begin(), fixMod.end());
      fixMod.resize(it_2 - fixMod.begin());
    }

    // append unassignedPeptideIdentifications
    unassigned_peptide_identifications_.insert(unassigned_peptide_identifications_.end(),
                                               rhs.unassigned_peptide_identifications_.begin(),
                                               rhs.unassigned_peptide_identifications_.end());

    // append consensusElements to consensusElementList:
    this->insert(this->end(), rhs.begin(), rhs.end());

    // todo: check for double entries
    // features, unassignedpeptides, proteins...

    // consistency
    try
    {
      UniqueIdIndexer<ConsensusMap>::updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition /*&e*/) // assign new UID's for conflicting entries
    {
      Size replaced_uids =  UniqueIdIndexer<ConsensusMap>::resolveUniqueIdConflicts();
      LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
    }

    return *this;
  }

  void ConsensusMap::clear(bool clear_meta_data)
  {
    Base::clear();

    if (clear_meta_data)
    {
      clearMetaInfo();
      clearRanges();
      // no "clear" method for DocumentIdentifier available
      this->DocumentIdentifier::operator=(DocumentIdentifier());
      clearUniqueId();
      file_description_.clear();
      experiment_type_.clear();
      protein_identifications_.clear();
      unassigned_peptide_identifications_.clear();
      data_processing_.clear();
    }
  }

  const ConsensusMap::FileDescriptions& ConsensusMap::getFileDescriptions() const
  {
    return file_description_;
  }

  ConsensusMap::FileDescriptions& ConsensusMap::getFileDescriptions()
  {
    return file_description_;
  }

  void ConsensusMap::setFileDescriptions(const ConsensusMap::FileDescriptions& file_description)
  {
    file_description_ = file_description;
  }

  const String& ConsensusMap::getExperimentType() const
  {
    return experiment_type_;
  }

  void ConsensusMap::setExperimentType(const String& experiment_type)
  {
    experiment_type_ = experiment_type;
  }

  void ConsensusMap::sortByIntensity(bool reverse)
  {
    if (reverse)
    {
      std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::IntensityLess()));
    }
    else
    {
      std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::IntensityLess());
    }
  }

  void ConsensusMap::sortByRT()
  {
    std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::RTLess());
  }

  void ConsensusMap::sortByMZ()
  {
    std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::MZLess());
  }

  void ConsensusMap::sortByPosition()
  {
    std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::PositionLess());
  }

  void ConsensusMap::sortByQuality(bool reverse)
  {
    if (reverse)
    {
      std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::QualityLess()));
    }
    else
    {
      std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::QualityLess());
    }
  }

  void ConsensusMap::sortBySize()
  {
    std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::SizeLess()));
  }

  void ConsensusMap::sortByMaps()
  {
    std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::MapsLess());
  }

  void ConsensusMap::sortPeptideIdentificationsByMapIndex()
  {
    // lambda predicate
    auto mapIndexLess = [] (const PeptideIdentification & a, const PeptideIdentification & b) -> bool
    {
      const bool has_a = a.metaValueExists("map_index");
      const bool has_b = b.metaValueExists("map_index");

      // moves IDs without meta value to end
      if (has_a && !has_b) { return true; }
      if (!has_a && has_b) { return false; }

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
        vector<PeptideIdentification> & pids = c.getPeptideIdentifications();
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
    Base::swap(from);

    // swap DocumentIdentifier
    DocumentIdentifier::swap(from);

    // swap unique id
    UniqueIdInterface::swap(from);

    // swap unique id index
    UniqueIdIndexer<ConsensusMap>::swap(from);

    // swap the remaining members
    std::swap(file_description_, from.file_description_);
    experiment_type_.swap(from.experiment_type_);
    protein_identifications_.swap(from.protein_identifications_);
    unassigned_peptide_identifications_.swap(from.unassigned_peptide_identifications_);
    data_processing_.swap(from.data_processing_);
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
    if (!s.empty())
    {
      this->setMetaValue("spectra_data", DataValue(s));
    }
  }

  /// get the file path to the first MS run
  void ConsensusMap::getPrimaryMSRunPath(StringList& toFill) const
  {
    if (this->metaValueExists("spectra_data"))
    {
      toFill = this->getMetaValue("spectra_data");
    }
  }

  /// Equality operator
  bool ConsensusMap::operator==(const ConsensusMap& rhs) const
  {
    return std::operator==(*this, rhs) &&
           MetaInfoInterface::operator==(rhs) &&
           RangeManagerType::operator==(rhs) &&
           DocumentIdentifier::operator==(rhs) &&
           UniqueIdInterface::operator==(rhs) &&
           file_description_ == rhs.file_description_ &&
           experiment_type_ == rhs.experiment_type_ &&
           protein_identifications_ == rhs.protein_identifications_ &&
           unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ &&
           data_processing_ == rhs.data_processing_;
  }

  /// Equality operator
  bool ConsensusMap::operator!=(const ConsensusMap& rhs) const
  {
    return !(operator==(rhs));
  }

  std::ostream& operator<<(std::ostream& os, const ConsensusMap& cons_map)
  {
    for (ConsensusMap::FileDescriptions::const_iterator it = cons_map.getFileDescriptions().begin(); it != cons_map.getFileDescriptions().end(); ++it)
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
    updateRanges_(begin(), end());

    // enlarge the range by the internal points of each feature
    for (Size i = 0; i < size(); ++i)
    {
      for (ConsensusFeature::HandleSetType::const_iterator it = operator[](i).begin(); it != operator[](i).end(); ++it)
      {
        double rt = it->getRT();
        double mz = it->getMZ();
        double intensity = it->getIntensity();

        // update RT
        if (rt < pos_range_.minPosition()[Peak2D::RT])
        {
          pos_range_.setMinX(rt);
        }
        if (rt > pos_range_.maxPosition()[Peak2D::RT])
        {
          pos_range_.setMaxX(rt);
        }
        // update m/z
        if (mz < pos_range_.minPosition()[Peak2D::MZ])
        {
          pos_range_.setMinY(mz);
        }
        if (mz > pos_range_.maxPosition()[Peak2D::MZ])
        {
          pos_range_.setMaxY(mz);
        }
        // update intensity
        if (intensity <  int_range_.minX())
        {
          int_range_.setMinX(intensity);
        }
        if (intensity > int_range_.maxX())
        {
          int_range_.setMaxX(intensity);
        }
      }
    }
  }

  bool ConsensusMap::isMapConsistent(Logger::LogStream* stream) const
  {
    Size stats_wrongMID(0); // invalid map ID references by a feature handle
    Map<Size, Size> wrong_ID_count; // which IDs were given which are not valid

    // check file descriptions
    std::set<String> maps;
    String all_maps; // for output later
    for (FileDescriptions::const_iterator it = file_description_.begin(); it != file_description_.end(); ++it)
    {
      String s = String("  file: ") + it->second.filename + " label: " + it->second.label;
      maps.insert(s);
      all_maps += s;
    }

    if (maps.size() != file_description_.size())
    {
      if (stream != nullptr)
      {
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
        if (file_description_.find(it->getMapIndex()) == file_description_.end())
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
        *stream << "ConsensusMap contains " << stats_wrongMID << " invalid references to maps:\n";
        for (Map<Size, Size>::ConstIterator it = wrong_ID_count.begin(); it != wrong_ID_count.end(); ++it)
        {
          *stream << "  wrong id=" << it->first << " (occurred " << it->second << "x)\n";
        }
        *stream << std::endl;
      }
      return false;
    }

    return true;
  }

} // namespace OpenMS
