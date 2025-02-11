// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

namespace OpenMS
{
  class AASequence;

  /**
    @brief This class implements a pair finding algorithm for consensus features.

    It offers a method to determine pairs across two consensus maps. The corresponding consensus
    features must be aligned, but may have small position deviations.

    The distance measure is implemented in class @ref FeatureDistance - see there for details.

    <B> Additional criteria for pairing </B>

    Depending on parameter @p use_identifications, peptide identifications annotated to the
    features may have to be compatible (i.e. no annotation or the same annotation) for a pairing
    to occur.

    Stability criterion: The distance to the nearest neighbor must be smaller than the distance
    to the second-nearest neighbor by a certain factor, see parameter @p second_nearest_gap.
    There is a non-trivial relation between this parameter and the maximum allowed difference
    (in RT or m/z) of the distance measure: If @p second_nearest_gap is greater than one,
    lowering @p max_difference may in fact lead to more - rather than fewer - pairings, because
    it increases the distance difference between the nearest and the second-nearest neighbor, so
    that the constraint imposed by @p second_nearest_gap may be fulfilled more often.

    <B> Quality calculation </B>

    The quality of a pairing is computed from the distance between the paired elements (nearest
    neighbors) and the distances to the second-nearest neighbors of both elements, according to
    the formula:

    @f[
    q_{i,j} = \big( 1 - d_{i,j} \big) \cdot
    \big( 1 - \frac{g \cdot d_{i,j}}{d_{2,i}} \big) \cdot
    \big( 1 - \frac{g \cdot d_{i,j}}{d_{2,j}} \big) \cdot
    @f]

    @f$ q_{i,j} @f$ is the quality of the pairing of elements @em i and @em j, @f$ d_{i,j} @f$ is
    the distance between the two, @f$ d_{2,i} @f$ and @f$d_{2,j} @f$ are the distances to the
    second-nearest neighbors of @em i and @em j, respectively, and @em g is the factor defined by
    parameter @p second_nearest_gap.

    Note that by the definition of the distance measure, @f$ 0 \leq d_{i,j} \leq 1 @f$ if @em i and
    @em j are to form a pair. The criteria for pairing further require that
    @f$ g \cdot d_{i,j} \leq d_{2,i} @f$ and @f$ g \cdot d_{i,j} \leq d_{2,j} @f$. This ensures that
    the resulting quality is always between one (best) and zero (worst).

    For the final quality @em q of the consensus feature produced by merging two paired elements
    (@em i and @em j), the existing quality values of the two elements are taken into account. The
    final quality is a weighted average of the existing qualities (@f$ q_i @f$ and @f$ q_j @f$) and
    the quality of the pairing (@f$ q_{i,j} @f$, see above):

    @f[
    q = \frac{q_{i,j} + (s_i - 1) \cdot q_i + (s_j - 1) \cdot q_j}{s_i + s_j - 1}
    @f]

    The weighting factors @f$ s_i @f$ and @f$ s_j @f$ are the sizes (i.e. numbers of subelements)
    of the two consensus features @em i and @em j. That way, it is possible to link several feature
    maps to a growing consensus map in a stepwise fashion (as done by @ref
    FeatureGroupingAlgorithmUnlabeled), and in the end obtain quality values that incorporate the
    qualities of all pairings that occurred during the generation of a consensus feature. Note that
    "missing" elements (if a consensus feature does not contain sub-features from all input maps)
    are not punished in this definition of quality.

    @htmlinclude OpenMS_StablePairFinder.parameters

    @ingroup FeatureGrouping
   */
  class OPENMS_DLLAPI StablePairFinder :
    public BaseGroupFinder
  {
public:

    ///Base class
    typedef BaseGroupFinder Base;

    /// Constructor
    StablePairFinder();

    /// Destructor
    ~StablePairFinder() override
    {
    }

    /**
      @brief Run the algorithm

      @note Exactly two @em input maps must be provided.

      @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<ConsensusMap>& input_maps,
             ConsensusMap& result_map) override;

protected:

    ///@name Internal helper classes and enums
    //@{
    enum
    {
      RT = Peak2D::RT,
      MZ = Peak2D::MZ
    };
    //@}

    //docu in base class
    void updateMembers_() override;

    /**
      @brief Checks if the peptide IDs of two features are compatible.

      A feature without identification is always compatible. Otherwise,
      two features are compatible if the best peptide hits of their
      identifications have the same sequences.
    */
    bool compatibleIDs_(const ConsensusFeature& feat1,
                        const ConsensusFeature& feat2) const;

    /// The distance to the second nearest neighbors must be by this factor larger than the distance to the matched element itself.
    double second_nearest_gap_;

    /// Only match if peptide IDs are compatible?
    bool use_IDs_;

    /**
      @brief Returns the highest scoring peptide hit in the the given peptide identification.

      @param peptideIdentification The peptideIdentification to scan.
    */
    const AASequence& getBestHitSequence_(const PeptideIdentification& peptideIdentification) const;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_STABLEPAIRFINDER_H

/*

 gnuplot history - how the plot was created - please do not delete this receipt

 f(x,intercept,exponent)=1/(1+(abs(x)*intercept)**exponent)
 set terminal postscript enhanced color
 set output "choosingstablepairfinderparams.ps"
 set size ratio .3
 plot [-3:3] [0:1] f(x,1,1), f(x,2,1), f(x,1,2), f(x,2,2)

 */
