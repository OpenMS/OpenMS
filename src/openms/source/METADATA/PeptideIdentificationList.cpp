// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  const String& PeptideIdentificationList::getListScoreType() const
  {
    return list_score_type_;
  }

  void PeptideIdentificationList::setListScoreType(const String& type)
  {
    list_score_type_ = type;
  }

  bool PeptideIdentificationList::isListHigherScoreBetter() const
  {
    return list_higher_score_better_;
  }

  void PeptideIdentificationList::setListHigherScoreBetter(bool value)
  {
    list_higher_score_better_ = value;
  }

  bool PeptideIdentificationList::usesListLevelScores() const
  {
    return use_list_level_scores_;
  }

  void PeptideIdentificationList::enableListLevelScores(bool enable)
  {
    use_list_level_scores_ = enable;
  }

  bool PeptideIdentificationList::hasConsistentScoreMetadata() const
  {
    if (empty()) return true;
    
    const String& ref_type = front().getScoreType();
    bool ref_orientation = front().isHigherScoreBetter();
    
    for (const auto& id : *this)
    {
      if (id.getScoreType() != ref_type || 
          id.isHigherScoreBetter() != ref_orientation)
      {
        return false;
      }
    }
    return true;
  }

  void PeptideIdentificationList::validateScoreConsistency() const
  {
    if (empty()) return;
    
    const String& ref_type = front().getScoreType();
    bool ref_orientation = front().isHigherScoreBetter();
    
    for (size_t i = 0; i < size(); ++i)
    {
      const auto& id = (*this)[i];
      if (id.getScoreType() != ref_type)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Inconsistent score types in PeptideIdentificationList at index " + String(i) + ": expected '" + ref_type + "'",
          id.getScoreType());
      }
      if (id.isHigherScoreBetter() != ref_orientation)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Inconsistent score orientations in PeptideIdentificationList at index " + String(i) + ": expected " + 
          String(ref_orientation ? "higher_better" : "lower_better"),
          String(id.isHigherScoreBetter() ? "higher_better" : "lower_better"));
      }
    }
  }

  void PeptideIdentificationList::migrateFromIndividualScores()
  {
    if (empty()) 
    {
      use_list_level_scores_ = true;
      return;
    }
    
    validateScoreConsistency();
    
    setListScoreType(front().getScoreType());
    setListHigherScoreBetter(front().isHigherScoreBetter());
    enableListLevelScores(true);
  }

  void PeptideIdentificationList::syncToIndividualScores()
  {
    if (!use_list_level_scores_) return;
    
    for (auto& id : *this)
    {
      id.setScoreType(list_score_type_);
      id.setHigherScoreBetter(list_higher_score_better_);
    }
  }

  void PeptideIdentificationList::tryMigrateAfterLoad(bool force_migration)
  {
    if (empty()) return;
    
    if (force_migration || hasConsistentScoreMetadata())
    {
      migrateFromIndividualScores();
      // Log the migration for debugging
      OPENMS_LOG_DEBUG << "Migrated PeptideIdentificationList to list-level score metadata: " 
                       << "score_type=" << list_score_type_ 
                       << ", higher_better=" << list_higher_score_better_ << std::endl;
    }
  }

  void PeptideIdentificationList::prepareForSave()
  {
    if (use_list_level_scores_)
    {
      syncToIndividualScores();
      // Log the sync for debugging
      OPENMS_LOG_DEBUG << "Synchronized list-level metadata to individual PeptideIdentification objects for file saving" << std::endl;
    }
  }
}
