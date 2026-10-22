// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Copies of the editable records of identification data, for an editor (e.g. the MetaDataBrowser); commit()
    writes the changed ones back at once

    The copies are the settings of runs, the database sequences of runs, the observations of identifications, the data
    of matches and the inference results. An editor changes them through the references these functions return (which
    stay valid as long as this object lives); commit() then replaces the records that differ from their originals.
    Records keep their IDs and runs their UUIDs, so links to them (e.g. of the features of a map) stay valid.

    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI IdentificationDataEdits
  {
  public:
    /// Edits of @p data, which must outlive this object
    explicit IdentificationDataEdits(IdentificationData& data);

    /// The identification data that commit() changes
    const IdentificationData& getData() const
    {
      return data_;
    }

    /// The copy of the settings of @p run (a run of the data)
    IdentificationData::RunSettings& settings(const IdentificationData::Run& run);

    /// The copy of the database sequences of @p run (a run of the data with database sequences)
    std::vector<IdentificationData::DatabaseSequence>& databaseSequences(const IdentificationData::Run& run);

    /// The copy of the observation of identification @p query of @p run
    IdentificationData::Observation& observation(const IdentificationData::Run& run, const IdentificationData::Identification& query);

    /// The copy of the data of match @p match of @p run (only payload edits such as metadata, see Run::replaceMatch())
    IdentificationData::MatchData& match(const IdentificationData::Run& run, const IdentificationData::Match& match);

    /// The copies of the inference results (in their order; add or remove none)
    std::vector<IdentificationData::InferenceResult>& inferenceResults()
    {
      return results_;
    }

    /**
      @brief Write the changed records to the data, all or none

      @return Whether anything changed
      @throw Exception::BaseException if a change is not valid (e.g. run settings that list 'spectra_data'); the data is
      unchanged then
    */
    bool commit();

  private:
    IdentificationData& data_;
    std::map<std::string, IdentificationData::RunSettings> settings_;
    std::map<std::string, std::vector<IdentificationData::DatabaseSequence>> sequences_;
    std::map<std::pair<std::string, UInt64>, IdentificationData::Observation> observations_;
    std::map<std::pair<std::string, UInt64>, IdentificationData::MatchData> matches_;
    std::vector<IdentificationData::InferenceResult> results_;
  };
} // namespace OpenMS
