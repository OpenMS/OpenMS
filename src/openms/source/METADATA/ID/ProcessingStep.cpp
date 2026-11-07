// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/ProcessingStep.h>

#include <tuple>

namespace OpenMS::IdentificationDataInternal
{
  bool ProcessingStep::operator<(const ProcessingStep& other) const
  {
    return (std::tie(date_time, software_ref, input_file_refs, actions) <
            std::tie(other.date_time, other.software_ref,
                     other.input_file_refs, other.actions));
  }

  bool ProcessingStep::operator==(const ProcessingStep& other) const
  {
    return (std::tie(software_ref, input_file_refs, date_time, actions) ==
            std::tie(other.software_ref, other.input_file_refs,
                     other.date_time, other.actions));
  }
} // namespace OpenMS::IdentificationDataInternal
