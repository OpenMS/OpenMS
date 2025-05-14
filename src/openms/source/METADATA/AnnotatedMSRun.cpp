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
  /**
   * @brief Returns a reference to the vector of peptide identifications.
   *
   * Allows direct access to the stored peptide identifications for modification or inspection.
   *
   * @return Reference to the internal vector of PeptideIdentification objects.
   */
  std::vector<PeptideIdentification>& AnnotatedMSRun::getPeptideIdentifications()
  {
    return peptide_ids_;
  }

  /**
   * @brief Returns a const reference to the vector of peptide identifications associated with this run.
   *
   * @return Const reference to the internal vector of PeptideIdentification objects.
   */
  const std::vector<PeptideIdentification>& AnnotatedMSRun::getPeptideIdentifications() const
  {
    return peptide_ids_;
  }

  /**
   * @brief Sets the peptide identifications for the run by copying from the provided vector.
   *
   * Replaces the current peptide identifications with a copy of the given vector.
   *
   * @param ids Vector of peptide identifications to assign.
   */
  void AnnotatedMSRun::setPeptideIdentifications(const std::vector<PeptideIdentification>& ids)
  {
    peptide_ids_ = ids;
  }

  /**
   * @brief Sets the peptide identifications by moving the provided vector.
   *
   * Transfers ownership of the given peptide identifications to the internal storage, enabling efficient data handling without copying.
   *
   * @param ids Rvalue reference to a vector of peptide identifications to be moved.
   */
  void AnnotatedMSRun::setPeptideIdentifications(std::vector<PeptideIdentification>&& ids)
  {
    peptide_ids_ = std::move(ids);
  }

  /**
   * @brief Returns a reference to the internal mass spectrometry experiment data.
   *
   * @return Reference to the stored MSExperiment object.
   */
  MSExperiment& AnnotatedMSRun::getMSExperiment()
  {
    return data;
  }

  /**
   * @brief Returns a const reference to the internal mass spectrometry experiment data.
   *
   * @return Const reference to the MSExperiment object managed by this run.
   */
  const MSExperiment& AnnotatedMSRun::getMSExperiment() const
  {
    return data;
  }

  /**
   * @brief Sets the internal mass spectrometry experiment data by moving from the provided experiment.
   *
   * Efficiently replaces the current MSExperiment with the given rvalue, transferring ownership of its contents.
   *
   * @param experiment Rvalue reference to the MSExperiment to be moved.
   */
  void AnnotatedMSRun::setMSExperiment(MSExperiment&& experiment)
  {
    data = std::move(experiment);
  }

  /**
   * @brief Sets the internal mass spectrometry experiment data by copying.
   *
   * Replaces the current MSExperiment with a copy of the provided experiment.
   *
   * @param experiment The MSExperiment object to copy.
   */
  void AnnotatedMSRun::setMSExperiment(const MSExperiment& experiment)
  {
    data = experiment;
  }
  
  /**
   * @brief Ensures the number of spectra matches the number of peptide identifications.
   *
   * Throws an Exception::InvalidValue if the internal MSExperiment and peptide identification counts differ, indicating an internal data inconsistency.
   *
   * @param function_name Name of the calling function for error reporting.
   */
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

