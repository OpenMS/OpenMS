// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{

  /**
    @brief A superimposer that uses a voting scheme, also known as pose clustering,
    to find a good affine transformation.

    This algorithm works on two consensus maps.  It computes an affine
    transformation that maps the elements of second map as near as possible
    to the elements in the first map.

    The voting scheme hashes affine transformations between pairs of features in
    map one and pairs of features in map two.  Each such pair of pairs defines a
    (potential) "pose" of the second map relative to the first.
    Then it finds a cluster in the parameter space of these poses.
    The affine transformation is then computed from this
    cluster of potential poses, hence the name pose clustering.

    @sa PoseClusteringShiftSuperimposer

    @htmlinclude OpenMS_PoseClusteringAffineSuperimposer.parameters

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI PoseClusteringAffineSuperimposer :
    public BaseSuperimposer
  {
public:

    /// Default ctor
    PoseClusteringAffineSuperimposer();

    /// Destructor
    
    ~PoseClusteringAffineSuperimposer() override
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

      @param map_model The model map (first input map)
      @param map_scene The scene map (second input map)
      @param transformation The output affine transformation (linear model transforming the scene map onto the model map)

      @exception IllegalArgument is thrown if the input maps are invalid.
    */
    void run(const ConsensusMap & map_model, const ConsensusMap & map_scene, TransformationDescription & transformation) override;

    /// Perform alignment on vector of 1D peaks
    virtual void run(const std::vector<Peak2D> & map_model, const std::vector<Peak2D> & map_scene, TransformationDescription & transformation);

  };
} // namespace OpenMS

