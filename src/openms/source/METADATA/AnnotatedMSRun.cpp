// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt, Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/METADATA/AnnotatedMSRun.h>


namespace OpenMS
{
  PeptideIdentificationList& AnnotatedMSRun::getPeptideIdentifications()
  {
    return peptide_ids_;
  }

  const PeptideIdentificationList& AnnotatedMSRun::getPeptideIdentifications() const
  {
    return peptide_ids_;
  }

  void AnnotatedMSRun::setPeptideIdentifications(const PeptideIdentificationList& ids)
  {
    peptide_ids_ = ids;
  }

  void AnnotatedMSRun::setPeptideIdentifications(PeptideIdentificationList&& ids)
  {
    peptide_ids_ = std::move(ids);
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
  
  void AnnotatedMSRun::checkPeptideIdSize_(const char* function_name) const
  {
    if (data.getSpectra().size() != peptide_ids_.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__,
        function_name, // Use the provided function name
        "Internal inconsistency: Number of spectra and peptide identifications do not match.",
        String(data.getSpectra().size()) + " vs " + String(peptide_ids_.size()));
    }
  }
}

