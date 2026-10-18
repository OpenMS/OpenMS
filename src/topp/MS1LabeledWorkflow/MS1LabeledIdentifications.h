// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <optional>
#include <vector>

namespace OpenMS
{
class MultiplexResolverAlgorithm;

/**
  @brief The identifications of the multiplets of one run, between their mapping and the resolver.

  The multiplets link their identifications as identification data: every assignment of an identification to a multiplet
  is an identification of its own, with the channel it matched as meta value @c map_index (see MS1LabeledFAIMS::annotate()).
*/
class MS1LabeledIdentifications
{
public:
  /// An identification to move: the multiplet that links it afterwards (none: unassigned) and its channel (meta value
  /// @c map_index; none: removed)
  struct Move
  {
    IdentificationData::QueryReference identification;
    std::optional<Size> multiplet;
    std::optional<Size> map_index;
  };

  /// The tolerances of the identification mapping (the workflow's section 'id_mapping')
  struct Tolerances
  {
    double rt;  ///< RT tolerance (seconds) beyond the elution region of a channel
    double mz;  ///< m/z tolerance
    bool ppm;   ///< @c mz in ppm (otherwise in Th)
  };

  /**
    @brief Move identifications of @p multiplets to other multiplets or to the unassigned ones.

    Every identification of @p moves is replaced by a copy with its matches and their scores, in the order of @p moves and
    after all other identifications, so that it comes last among the identifications of its multiplet resp. among the
    unassigned ones (like a peptide identification appended to a list). The links to the original are removed.

    @throw Exception::InvalidValue if an identification is unknown or to be moved twice (nothing changes)
  */
  static void moveIdentifications(ConsensusMap& multiplets, const std::vector<Move>& moves);

  /**
    @brief Attach identifications of channels the detection did not find to the multiplet of their peptide.

    An identification maps onto a multiplet through its own channel. When that channel was too weak to be detected (e.g. the
    minor channel of a mixture, which DDA fragments all the same), the multiplet can still be there, detected in the other
    channels of the peptide. So an identification that stayed unassigned is looked up at the m/z its peptide has in the other
    channels (MultiplexResolverAlgorithm::expectedChannelMZs()): if a partial multiplet of the same charge and CV has a channel
    there (one of @p channels, see the workflow) that elutes at the identification (within @p tolerances), the identification
    is attached to it, and its own channel is added to the multiplet as a feature of intensity 0 for the re-quantification to
    measure. The resolver then checks the result like any other identified multiplet.

    @param[in,out] multiplets Multiplets with their identifications
    @param[in] channels One feature per channel of a multiplet, with the unique id of its feature handle and its elution region as convex hulls
    @param[in] resolver Gives the m/z of a peptide in every channel
    @param[in] multiplicity Number of channels: a multiplet with that many detected channels is complete
    @param[in] tolerances Tolerances of the identification mapping
    @return The number of identifications attached
  */
  static Size attachToUndetectedChannels(ConsensusMap& multiplets, const FeatureMap& channels, const MultiplexResolverAlgorithm& resolver,
                                         Size multiplicity, const Tolerances& tolerances);

  /**
    @brief Attach every spectrum match to the closest channel, and to the closest multiplet only.

    The mapping attaches an identification to every multiplet within its tolerances and records the channel it matched as its
    @c map_index. The resolver reads that index as the channel the spectrum belongs to, so it is set to the channel closest to
    the identification's m/z here: the theoretical m/z of the identified peptidoform, whose labels name its channel (the
    precursor m/z can be off by isotopes and then lie closer to another channel, e.g. Lys4 at charge 3). One spectrum is
    evidence for one multiplet, too (and a QPX PSM may reference one feature): it stays with the multiplet whose closest
    channel is closest to the identification in m/z, then in RT. This runs before the best-scoring match is chosen per
    multiplet, so that a multiplet does not first take a neighbour's better-scoring spectrum and then lose it, ending up
    without its own. Distinct spectra of the same peptide on distinct multiplets are not touched: both multiplets are
    quantified and their channel values add up per peptide.

    @return The number of spectrum matches (spectrum, peptidoform, charge) that were on several multiplets
  */
  static Size keepMultiAssignedOnClosest(ConsensusMap& multiplets);

  /**
    @brief Of several identifications on one multiplet, keep those whose labels agree with its mass shifts.

    The best-scoring identification is the one a multiplet is quantified under, and if its labels contradict the detected mass
    shifts, the multiplet is a quant/ID conflict - even when another identification on it fits (e.g. the best one comes from
    a co-isolated peptide). So the identifications that fit (MultiplexResolverAlgorithm::checkIdentifications()) are
    preferred: those that do not are moved to the unassigned ones (they remain PSMs for protein inference), without channel.
    A multiplet that none of its identifications fits stays a conflict.

    @return The number of identifications moved to the unassigned ones
  */
  static Size keepConsistentIdentifications(ConsensusMap& multiplets, const MultiplexResolverAlgorithm& resolver);
};
} // namespace OpenMS
