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
  PeptideIdentification& AnnotatedMSRawData::getPeptideIdentification(size_t index)
  {
    return peptide_ids[index];
  }

  const PeptideIdentification& AnnotatedMSRawData::getPeptideIdentification(size_t index) const
  {
    return peptide_ids[index];
  }

  std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications()
  {
    return peptide_ids;
  }

  const std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications() const
  {
    return peptide_ids;
  }

  void AnnotatedMSRawData::setPeptideIdentification(PeptideIdentification&& id, size_t index)
  {
    peptide_ids[index] = std::move(id);
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
