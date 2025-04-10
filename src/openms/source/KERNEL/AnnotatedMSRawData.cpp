// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/AnnotatedMSRawData.h>

namespace OpenMS
{
  std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications()
  {
    return peptide_ids;
  }

  const std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications() const
  {
    return peptide_ids;
  }

  void AnnotatedMSRawData::setPeptideIdentifications(std::vector<PeptideIdentification>&& ids)
  {
    peptide_ids = std::move(ids);
  }

  MSExperiment& AnnotatedMSRawData::getMSExperiment()
  {
    return data;
  }

  const MSExperiment& AnnotatedMSRawData::getMSExperiment() const
  {
    return data;
  }
}
