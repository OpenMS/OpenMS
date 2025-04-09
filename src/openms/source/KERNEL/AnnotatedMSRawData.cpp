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
  std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications(size_t index)
  {
    return peptide_ids[index];
  }

  const std::vector<PeptideIdentification>& AnnotatedMSRawData::getPeptideIdentifications(size_t index) const
  {
    return peptide_ids[index];
  }

  std::vector<std::vector<PeptideIdentification>>& AnnotatedMSRawData::getAllPeptideIdentifications()
  {
    return peptide_ids;
  }

  const std::vector<std::vector<PeptideIdentification>>& AnnotatedMSRawData::getAllPeptideIdentifications() const
  {
    return peptide_ids;
  }

  void AnnotatedMSRawData::setPeptideIdentifications(std::vector<PeptideIdentification>&& ids, size_t index)
  {
    peptide_ids[index] = std::move(ids);
  }

  void AnnotatedMSRawData::setAllPeptideIdentifications(std::vector<std::vector<PeptideIdentification>>&& ids)
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
