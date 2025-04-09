// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <map>
#include <vector>
namespace OpenMS
{
  class FeatureMap;
  /**
   * @brief This class is a metric for the QualityControl TOPP Tool.
   *
   * This class counts the number of MissedCleavages per PeptideIdentification given a FeatureMap
   * and returns an agglomeration statistic (observed counts).
   * Additionally the PeptideHits in the FeatureMap are augmented with MetaInformation:
   *  - 'missed_cleavages'
   *  - 'FWHM' (from feature's 'FWHM' or 'model_FWHM')
   *  - 'mass' (experimental mass of peptide)
   */
  class OPENMS_DLLAPI MissedCleavages : public QCBase
  {
  private:
    typedef std::map<UInt32, UInt32> MapU32;
    /// collects number of missed cleavages from PeptideIdentification in a result map (missed cleavages: occurences)
    void get_missed_cleavages_from_peptide_identification_(const ProteaseDigestion& digestor, MapU32& result, const UInt32& max_mc, PeptideIdentification& pep_id);

  public:
    /// constructor
    MissedCleavages() = default;

    /// destructor
    virtual ~MissedCleavages() = default;

    /**
     * @brief Counts the number of missed cleavages per PeptideIdentification.
     *
     * The result is a key/value map: \#missed_cleavages --> counts
     * Additionally the first PeptideHit in each PeptideIdentification of the FeatureMap is annotated with metavalue 'missed_cleavages'.
     * The protease and digestion parameters are taken from the first ProteinIdentication (and SearchParameter therein) within the FeatureMap itself.
     *
     * @param fmap FeatureMap with Peptide and ProteinIdentifications
     */
    void compute(FeatureMap& fmap);
    void compute(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    /// returns the name of the metric
    const String& getName() const override;

    /// returns the result as maps of number of missed_cleavages to counts; one map for each call to compute(...)
    const std::vector<std::map<UInt32, UInt32>>& getResults() const;

    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for POSTFDRFEAT;
     */
    QCBase::Status requirements() const override;

  private:
    /// container that stores results
    std::vector<std::map<UInt32, UInt32>> mc_result_;
  };
} // namespace OpenMS
