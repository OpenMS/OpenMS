// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <set>
#include <map>

namespace OpenMS
{

  class MassExplainer;
  class FeatureMap;
  class ChargePair;

  /**
    @brief ILP-based decharging: pick the maximum-log-score subset of charge-pair edges that yields a consistent per-feature charge assignment.

    Wraps the integer-linear-programming step shared by @ref FeatureDeconvolution (peptides)
    and @ref MetaboliteFeatureDeconvolution (metabolites). Given a @ref FeatureMap and a
    pre-computed list of putative @ref ChargePair edges between its features (each edge
    captures one consistent charge-state hypothesis for a pair of features that are also
    adduct/mass-shift compatible), the solver returns the maximum-log-score subset of edges
    that respects per-feature consistency constraints (at most one realised
    charge-assignment per feature). The realised edges are marked active in @p pairs.

    The implementation first partitions the input edge set into connected components and
    solves one ILP per component (so a single big input doesn't translate into a single
    intractable LP) — this slicing is handled internally and is transparent to callers.

    LP-solver choice (GLPK, COIN-OR, or HiGHS — selected at build time via the
    @c LP_SOLVER CMake variable) is hidden behind @ref LPWrapper.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI ILPDCWrapper
  {

public:
    /// Container of putative charge-pair edges between features
    typedef std::vector<ChargePair> PairsType;
    /// Index type into @ref PairsType, used for slicing
    typedef PairsType::size_type PairsIndex;

    /// Default constructor; the class is stateless and only exposes compute()
    ILPDCWrapper();

    /// Destructor
    virtual ~ILPDCWrapper();

    /**
      @brief Solve the ILP on @p pairs to find the max-log-score subset of consistent charge-pair edges.

      Partitions the input edge graph into connected components, solves one ILP per
      component (per @ref LPWrapper), accumulates the objective values, and marks every
      realised edge active in @p pairs. The @ref FeatureMap is read-only here — features
      themselves are not modified; the caller is expected to apply the chosen charge
      assignments to @p fm separately.

      @param[in]     fm            Feature map providing per-feature mass / intensity for log-score evaluation.
      @param[in,out] pairs         Putative charge-pair edges to filter; realised edges are flagged active in place.
      @param[in]     verbose_level Verbosity for the LP solver (0 = silent; higher = more solver chatter).
      @return Sum of per-component objective values, or @c -1 if @p fm is empty (a warning is logged in that case).
    */
    double compute(const FeatureMap& fm, PairsType& pairs, Size verbose_level) const;

private:

    /// Solve one slice of the edge graph (the new clique-based solver) — see compute() for the externally visible contract
    double computeSlice_(const FeatureMap& fm,
                         PairsType& pairs,
                         const PairsIndex margin_left,
                         const PairsIndex margin_right,
                         const Size verbose_level) const;

    /// Legacy slice solver kept for comparison / regression testing — kept in lock-step with @ref computeSlice_
    double computeSliceOld_(const FeatureMap& fm,
                            PairsType& pairs,
                            const PairsIndex margin_left,
                            const PairsIndex margin_right,
                            const Size verbose_level) const;

    /// Log-score of one putative charge-pair edge given the underlying feature intensities and the recorded mass-shift probability
    double getLogScore_(const PairsType::value_type& pair, const FeatureMap& fm) const;

    /// Internal lookup: feature id (as @ref String) -> set of compatible charge-annotation variant indices encountered for that feature so far
    typedef std::map<String, std::set<Size> > FeatureType_;

    /// Record one more charge-annotation variant @p v for the feature whose rotated id is @p rota_l in @p f_set
    void updateFeatureVariant_(FeatureType_& f_set, const String& rota_l, const Size& v) const;



  }; // !class

} // !namespace

