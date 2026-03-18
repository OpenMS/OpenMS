// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

namespace OpenMS
{

  /**
   @brief  A superimposer that uses a voting scheme, also known as pose clustering,
   to find a good shift transformation.

   This algorithm works on two consensus maps.  It computes an shift
   transformation that maps the elements of second map as near as possible
   to the elements in the first map.

   The voting scheme hashes shift transformations between features in
   map one and features in map two.  Each such pair defines a
   (potential) "pose" of the second map relative to the first.
   Then it finds a cluster in the parameter space of these poses.
   The shift transformation is then computed from this
   cluster of potential poses, hence the name pose clustering.

   @sa PoseClusteringAffineSuperimposer

   @htmlinclude OpenMS_PoseClusteringShiftSuperimposer.parameters

   @ingroup MapAlignment
   */
  class OPENMS_DLLAPI PoseClusteringShiftSuperimposer :
    public BaseSuperimposer
  {
public:

    /// Default ctor
    PoseClusteringShiftSuperimposer();

    /// Destructor
    ~PoseClusteringShiftSuperimposer() override
    {}

    /**
      @brief Estimates the transformation and fills the given mapping function. (Has a precondition!)

      @note Exactly two input maps must be given.

      @pre For performance reasons, we trust that (the equivalent of:)
      <code>
        maps[0].updateRanges();
        maps[1].updateRanges();
      </code>
      has been done <i>before</i> calling this.  You have been warned!

      @exception IllegalArgument is thrown if the input maps are invalid.
    */
    void run(const ConsensusMap & map_model, const ConsensusMap & map_scene, TransformationDescription & transformation) override;

  };
} // namespace OpenMS

