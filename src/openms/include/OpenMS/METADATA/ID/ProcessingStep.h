// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ID/ProcessingSoftware.h>
#include <OpenMS/METADATA/ID/InputFile.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct ProcessingStep: public MetaInfoInterface
    {
      ProcessingSoftwareRef software_ref;

      std::vector<InputFileRef> input_file_refs;

      DateTime date_time;

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      explicit ProcessingStep(
        ProcessingSoftwareRef software_ref,
        const std::vector<InputFileRef>& input_file_refs = std::vector<InputFileRef>(),
        const DateTime& date_time = DateTime::now(),
        const std::set<DataProcessing::ProcessingAction>& actions = std::set<DataProcessing::ProcessingAction>())
        :
        software_ref(software_ref), input_file_refs(input_file_refs),
        date_time(date_time), actions(actions)
      {
      }

      ProcessingStep(const ProcessingStep& other) = default;

      // order by date/time first, don't compare meta data (?):
      bool operator<(const ProcessingStep& other) const
      {
        return (std::tie(date_time, software_ref, input_file_refs, actions) <
                std::tie(other.date_time, other.software_ref,
                         other.input_file_refs, other.actions));
      }

      // don't compare meta data (?):
      bool operator==(const ProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, date_time, actions) ==
                std::tie(other.software_ref, other.input_file_refs,
                         other.date_time, other.actions));
      }
    };

    typedef std::set<ProcessingStep> ProcessingSteps;
    typedef IteratorWrapper<ProcessingSteps::iterator> ProcessingStepRef;

  }
}
