// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <functional>
#include <string>
#include <vector>

namespace OpenMS
{
  class PeptideHit;

  /**
    @brief FDR within classes of PSMs (e.g. modified vs. unmodified, isoform-specific vs. canonical)

    One q-value over all PSMs lets a small class pass at a much higher error rate than the large class: its decoys
    are diluted among the large class's targets. ClassFDR gives every target PSM a q-value within its class; a PSM is
    then accepted at a cutoff only if both its global and its class q-value pass it (see raiseToClassQValues()).

    Within a class, each target gets an empirical p-value against all decoys (the pooled null), and Storey q-values
    (lambda = 0.5) are computed over the class's targets. This borrows the power of the large decoy set and is
    conservative for a class whose own decoys score no higher than the pool. That is checked: a class with fewer than
    20 decoys, or whose decoys fall in the pool's top 5% significantly more often than 5% (one-sided z > 2.33), uses
    target-decoy competition inside the class instead, q = (decoys + 1) / targets.

    Ported from bigbio/andes #164.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ClassFDR
  {
  public:
    /// how the q-values of a class were estimated
    enum class Method
    {
      POOLED_NULL, ///< pooled-decoy p-values, Storey q-values within the class
      CLASS_TDC    ///< target-decoy competition within the class
    };

    /// per-class summary
    struct OPENMS_DLLAPI ClassSummary
    {
      std::string name;
      Size targets = 0;
      Size decoys = 0;
      Method method = Method::POOLED_NULL;
      double pi0 = 1.0;
    };

    /// meta value holding the class q-value of a hit
    static const std::string CLASS_Q_VALUE;

    /**
      @brief Annotate every target hit with its class q-value (meta value CLASS_Q_VALUE)

      Must run before decoys are removed. Hits are ranked by the main score of their identification (its orientation
      from PeptideIdentification::isHigherScoreBetter()). Decoys are hits with target_decoy "decoy". With a single class,
      nothing is annotated.

      @param ids Identifications
      @param class_of The class of a hit
      @param use_all_hits Use all hits, not only the best hit of each identification
      @return The classes, largest first
    */
    static std::vector<ClassSummary> annotateClassQValues(PeptideIdentificationList& ids,
                                                          const std::function<std::string(const PeptideHit&)>& class_of,
                                                          bool use_all_hits = false);

    /**
      @brief Raise each hit's q-value (main score, which must be a q-value) to its class q-value, where annotated
      @return Number of hits whose q-value was raised
    */
    static Size raiseToClassQValues(PeptideIdentificationList& ids);

    /// Class of a hit from the meta value @p name (its string form; "NA" if missing)
    static std::string classFromMetaValue(const PeptideHit& hit, const std::string& name);

    /// Class of a hit from its modifications: "unmodified", or the sorted names of its modifications joined by '+',
    /// ignoring those in @p ignore (e.g. fixed modifications)
    static std::string classFromModifications(const PeptideHit& hit, const std::vector<std::string>& ignore);

    /// Storey q-values (lambda = 0.5) of @p p in the same order; @p pi0 receives the null-fraction estimate
    static std::vector<double> storeyQValues(const std::vector<double>& p, double& pi0);
  };
} // namespace OpenMS
