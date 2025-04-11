// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/AnnotatedMSRun.h>

namespace OpenMS
{
  std::vector<PeptideIdentification>& AnnotatedMSRun::getPeptideIdentifications()
  {
    return peptide_ids;
  }

  const std::vector<PeptideIdentification>& AnnotatedMSRun::getPeptideIdentifications() const
  {
    return peptide_ids;
  }

  void AnnotatedMSRun::setPeptideIdentifications(const std::vector<PeptideIdentification>& ids)
  {
    peptide_ids = ids;
  }

  void AnnotatedMSRun::setPeptideIdentifications(std::vector<PeptideIdentification>&& ids)
  {
    peptide_ids = std::move(ids);
  }

  MSExperiment& AnnotatedMSRun::getMSExperiment()
  {
    return data;
  }

  const MSExperiment& AnnotatedMSRun::getMSExperiment() const
  {
    return data;
  }

  void AnnotatedMSRun::setMSExperiment(MSExperiment&& experiment)
  {
    data = std::move(experiment);
  }

  void AnnotatedMSRun::setMSExperiment(const MSExperiment& experiment)
  {
    data = experiment;
  }

}
