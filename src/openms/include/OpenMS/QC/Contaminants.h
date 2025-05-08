// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <unordered_set>


namespace OpenMS
{
  /**
   * @brief This class is a metric for the QualityControl TOPP tool.
   *
   * This class checks whether a peptide is a contaminant (given a protein DB) and adds that result as metavalue "is_contaminant"
   * to the first hit of each PeptideIdentification.
   */
  class OPENMS_DLLAPI Contaminants : public QCBase
  {
  public:
    /// structure for storing results
    struct ContaminantsSummary {
      ///(\#contaminants in assigned/ \#peptides in assigned)
      double assigned_contaminants_ratio;

      ///(\#contaminants in unassigned/ \#peptides in unassigned)
      double unassigned_contaminants_ratio;

      ///(\#all contaminants/ \#peptides in all)
      double all_contaminants_ratio;

      ///(intensity of contaminants in assigned/ intensity of peptides in assigned)
      double assigned_contaminants_intensity_ratio;

      ///(features without peptideidentification or with peptideidentifications but without hits; all features)
      std::pair<Int64, Int64> empty_features;
    };

    /// Constructor
    Contaminants() = default;

    /// Destructor
    virtual ~Contaminants() = default;

    /**
     * @brief Checks if the peptides are in the contaminant database.
     *
     * "is_contaminant" metavalue is added to the first hit of each PeptideIdentification of each feature
     * and to the first hit of all unsigned PeptideIdentifications.
     * The enzyme and number of missed cleavages used to digest the given protein DB is taken
     * from the ProteinIdentification[0].getSearchParameters() within the given FeatureMap.
     *
     * @param features Input FeatureMap with peptideidentifications of features
     * @param contaminants Vector of FASTAEntries that need to be digested to check whether a peptide is a contaminant or not
     * @exception Exception::MissingInformation if the contaminants database is empty
     * @exception Exception::MissingInformation if no enzyme is given
     * @exception Exception::MissingInformation if proteinidentification of FeatureMap is empty
     * @warning LOG_WARN if the FeatureMap is empty
     */
    void compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants);

    /// returns the name of the metric
    const String& getName() const override;

    /// returns results
    const std::vector<Contaminants::ContaminantsSummary>& getResults();

    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for POSTFDRFEAT and CONTAMINANTS
     */
    Status requirements() const override;

  private:
    /// name of the metric
    const String name_ = "Contaminants";

    /// container that stores results
    std::vector<Contaminants::ContaminantsSummary> results_;

    /// unordered set that contains the contaminant sequences
    std::unordered_set<String> digested_db_;

    /**
     * @brief
     * checks if the peptide is in the contaminant database
     * @param key String that will be the key for searching in the unordered set
     * @param pep_hit PeptideHit to store the result "is_contaminant = 0/1"
     * @param total counter of all checked peptides
     * @param cont counter of all checked peptides that are contaminants
     * @param sum_total intensity of all checked peptides
     * @param sum_cont intensity of all checked peptides that are contaminants
     * @param intensity intensity of current peptide
     */
    void compare_(const String& key, PeptideHit& pep_hit, Int64& total, Int64& cont, double& sum_total, double& sum_cont, double intensity);
  };
} // namespace OpenMS
