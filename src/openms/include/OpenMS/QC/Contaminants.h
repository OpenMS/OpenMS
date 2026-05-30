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
    @brief QualityControl metric: flag PSMs whose peptide sequences match
           a user-supplied contaminants FASTA (e.g. cRAP) after digestion.

    Each call to @ref compute marks the first @c PeptideHit of every
    @c PeptideIdentification in the supplied @c FeatureMap (both
    feature-attached and unassigned) with an @c "is_contaminant"
    meta value (@c 1 if its unmodified sequence matches an entry in
    the digested contaminants database, @c 0 otherwise) and appends a
    @ref ContaminantsSummary to the internal results list returned by
    @ref getResults.

    @note The digested contaminants database is cached on the first
          @ref compute call; subsequent calls reuse the same cache
          and ignore the @p contaminants argument's actual contents.
          The enzyme and missed-cleavage settings used for the initial
          digestion are also frozen at that point (taken from
          @c FeatureMap::getProteinIdentifications()[0].getSearchParameters()
          of the @c FeatureMap supplied on the first call).

    @ingroup Metadata
  */
  class OPENMS_DLLAPI Contaminants : public QCBase
  {
  public:
    /**
      @brief Result bundle returned per @ref compute call.

      All four ratios are unitless and lie in @c [0, 1] for non-empty
      denominators; when the corresponding denominator is @c 0 the
      ratio evaluates to @c NaN.
    */
    struct ContaminantsSummary {
      double assigned_contaminants_ratio;            ///< \#contaminants in feature-attached PSMs / \#feature-attached PSMs.
      double unassigned_contaminants_ratio;          ///< \#contaminants in unassigned PSMs / \#unassigned PSMs.
      double all_contaminants_ratio;                 ///< \#contaminants in all PSMs / \#all PSMs.
      double assigned_contaminants_intensity_ratio;  ///< Sum of feature intensities of contaminant feature-attached PSMs / sum of feature intensities of all feature-attached PSMs.
      std::pair<Int64, Int64> empty_features;        ///< (Number of features without a peptide hit, total number of features).
    };

    /// Default constructor.
    Contaminants() = default;

    /// Destructor.
    virtual ~Contaminants() = default;

    /**
      @brief Annotate the PSMs of @p features with @c "is_contaminant" and append a summary to @ref getResults.

      On the first call (when the internal digested-DB cache is empty),
      the FASTA entries in @p contaminants are digested using the
      enzyme and missed-cleavage count from
      @c features.getProteinIdentifications()[0].getSearchParameters()
      and stored in a hash set. Every subsequent call to @ref compute
      reuses that cache without consulting @p contaminants again.

      The first @c PeptideHit of every @c PeptideIdentification in
      @p features (both attached to features and in
      @c getUnassignedPeptideIdentifications()) is annotated with
      @c "is_contaminant" set to @c 0 or @c 1. The aggregated ratios
      and the empty-feature counters are then pushed onto the internal
      results list (see @ref getResults).

      @param[in,out] features    Source of the PSMs (annotated in place)
                                 and -- on the first call -- of the
                                 digestion enzyme used to build the
                                 cache.
      @param[in]     contaminants Contaminants FASTA. Only consulted on
                                  the first call; subsequent calls
                                  ignore it.

      @throws Exception::MissingInformation when @p contaminants is
                                            empty.
      @throws Exception::MissingInformation when the FeatureMap has no
                                            protein identification
                                            (only on the first call).
      @throws Exception::MissingInformation when the configured
                                            digestion enzyme is
                                            @c "unknown_enzyme" (only
                                            on the first call).

      @note When @p features is empty, a warning is logged and the
            method still runs (and may produce @c NaN ratios).
    */
    void compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants);

    /// Name of this QC metric (@c "Contaminants").
    const String& getName() const override;

    /// Per-call summaries appended by @ref compute, in call order.
    const std::vector<Contaminants::ContaminantsSummary>& getResults();

    /**
      @brief Input-data requirements of @ref compute.

      @return Status flags with @c POSTFDRFEAT and @c CONTAMINANTS set.
    */
    Status requirements() const override;

  private:
    /// Metric name returned by @ref getName.
    const String name_ = "Contaminants";

    /// Per-call summaries; @ref compute appends one entry per invocation.
    std::vector<Contaminants::ContaminantsSummary> results_;

    /// Cached digested contaminants database, filled on the first @ref compute call and reused thereafter.
    std::unordered_set<String> digested_db_;

    /**
      @brief Increment the contaminant counters and annotate one hit with @c "is_contaminant".

      Looks @p key up in @ref digested_db_ and updates the counters
      and @p pep_hit accordingly.

      @param[in]     key       Unmodified peptide sequence to look up.
      @param[in,out] pep_hit   Hit to annotate with @c "is_contaminant" set to @c 0 (not found) or @c 1 (found).
      @param[in,out] total     Total-PSMs counter; incremented by @c 1.
      @param[in,out] cont      Contaminant-PSMs counter; incremented only on a hit.
      @param[in,out] sum_total Total-intensity accumulator; incremented by @p intensity.
      @param[in,out] sum_cont  Contaminant-intensity accumulator; incremented by @p intensity only on a hit.
      @param[in]     intensity Intensity associated with this PSM.
    */
    void compare_(const String& key, PeptideHit& pep_hit, Int64& total, Int64& cont, double& sum_total, double& sum_cont, double intensity);
  };
} // namespace OpenMS
