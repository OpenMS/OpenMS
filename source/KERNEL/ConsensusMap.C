// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusMap.h>

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

  ConsensusMap::ConsensusMap(const ConsensusMap & source) :
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

  ConsensusMap & ConsensusMap::operator=(const ConsensusMap & source)
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

  ConsensusMap & ConsensusMap::operator+=(const ConsensusMap & rhs)
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

    // append proteinIdenficiation
    protein_identifications_.insert(protein_identifications_.end(),
                                    rhs.protein_identifications_.begin(),
                                    rhs.protein_identifications_.end());

    // ensure non-redundant modification parameter
    for (std::vector<ProteinIdentification>::iterator it = protein_identifications_.begin();
         it != protein_identifications_.end();
         ++it)
    {
      std::vector<String>::iterator it2;

      // remove redundant variable modifications
      std::vector<String> & varMod = const_cast<std::vector<String> &>(it->getSearchParameters().variable_modifications);
      sort(varMod.begin(), varMod.end());
      it2 = unique(varMod.begin(), varMod.end());
      varMod.resize(it2 - varMod.begin());

      // remove redundant fixed modifications
      std::vector<String> & fixMod = const_cast<std::vector<String> &>(it->getSearchParameters().fixed_modifications);
      sort(fixMod.begin(), fixMod.end());
      it2 = unique(fixMod.begin(), fixMod.end());
      fixMod.resize(it2 - fixMod.begin());
    }

    // append unassignedPeptideIdentifiactions
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

  const ConsensusMap::FileDescriptions & ConsensusMap::getFileDescriptions() const
  {
    return file_description_;
  }

  ConsensusMap::FileDescriptions & ConsensusMap::getFileDescriptions()
  {
    return file_description_;
  }

  const String & ConsensusMap::getExperimentType() const
  {
    return experiment_type_;
  }

  void ConsensusMap::setExperimentType(const String & experiment_type)
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

  void ConsensusMap::swap(ConsensusMap & from)
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
  const std::vector<ProteinIdentification> & ConsensusMap::getProteinIdentifications() const
  {
    return protein_identifications_;
  }

  /// mutable access to the protein identifications
  std::vector<ProteinIdentification> & ConsensusMap::getProteinIdentifications()
  {
    return protein_identifications_;
  }

  /// sets the protein identifications
  void ConsensusMap::setProteinIdentifications(const std::vector<ProteinIdentification> & protein_identifications)
  {
    protein_identifications_ = protein_identifications;
  }

  /// non-mutable access to the unassigned peptide identifications
  const std::vector<PeptideIdentification> & ConsensusMap::getUnassignedPeptideIdentifications() const
  {
    return unassigned_peptide_identifications_;
  }

  /// mutable access to the unassigned peptide identifications
  std::vector<PeptideIdentification> & ConsensusMap::getUnassignedPeptideIdentifications()
  {
    return unassigned_peptide_identifications_;
  }

  /// sets the unassigned peptide identifications
  void ConsensusMap::setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification> & unassigned_peptide_identifications)
  {
    unassigned_peptide_identifications_ = unassigned_peptide_identifications;
  }

  /// returns a const reference to the description of the applied data processing
  const std::vector<DataProcessing> & ConsensusMap::getDataProcessing() const
  {
    return data_processing_;
  }

  /// returns a mutable reference to the description of the applied data processing
  std::vector<DataProcessing> & ConsensusMap::getDataProcessing()
  {
    return data_processing_;
  }

  /// sets the description of the applied data processing
  void ConsensusMap::setDataProcessing(const std::vector<DataProcessing> & processing_method)
  {
    data_processing_ = processing_method;
  }

  /// Equality operator
  bool ConsensusMap::operator==(const ConsensusMap & rhs) const
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
  bool ConsensusMap::operator!=(const ConsensusMap & rhs) const
  {
    return !(operator==(rhs));
  }

  std::ostream & operator<<(std::ostream & os, const ConsensusMap & cons_map)
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

    //enlarge the range by the internal points of each feature
    for (Size i = 0; i < size(); ++i)
    {
      for (ConsensusFeature::HandleSetType::const_iterator it = operator[](i).begin(); it != operator[](i).end(); ++it)
      {
        DoubleReal rt = it->getRT();
        DoubleReal mz = it->getMZ();
        DoubleReal intensity = it->getIntensity();

        //update RT
        if (rt < pos_range_.minPosition()[Peak2D::RT])
        {
          pos_range_.setMinX(rt);
        }
        if (rt > pos_range_.maxPosition()[Peak2D::RT])
        {
          pos_range_.setMaxX(rt);
        }
        //update m/z
        if (mz < pos_range_.minPosition()[Peak2D::MZ])
        {
          pos_range_.setMinY(mz);
        }
        if (mz > pos_range_.maxPosition()[Peak2D::MZ])
        {
          pos_range_.setMaxY(mz);
        }
        //update intensity
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

} // namespace OpenMS
