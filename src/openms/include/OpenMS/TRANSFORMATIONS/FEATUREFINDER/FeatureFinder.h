// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**@brief The main feature finder class.

      - Stores the flags for (indices of) data points ("used", "unused")
      - The algorithm itself is a factory product (derived from FeatureFinderAlgorithm)
      - The main method is run(), which is a template so that we can deal with different types of input and output
      - The run() method takes five arguments: algorithm_name, input_map, output, parameters, seeds
      .

      @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI FeatureFinder :
    public ProgressLogger,
    public FeatureFinderDefs
  {

public:
    /// Default constructor.
    FeatureFinder();

    /// Destructor
    ~FeatureFinder() override;

    /**
        @brief Executes the FeatureFinder using the given algorithm

        There are several constraints for the @p input_map.  They are tested before
        the algorithm starts.  It must only contain MS 1 level scans and you
        have to call updateRanges() before passing it to this method.
The input map is sorted by RT & m/z if that's not the case.
Furthermore we throw an Exception if the data contains negative m/z values,
as this will disturb most algorithms.

        @param algorithm_name Name of the feature finding algorithm to use
        @param input_map Input peak map
        @param features Output feature map
        @param param Algorithm parameters
        @param seeds List of seeds to use

        Implemented in FeatureFinder_impl.h
    */
    void run(const String& algorithm_name, PeakMap& input_map, FeatureMap& features, const Param& param, const FeatureMap& seeds);

    /// Returns a non-mutable reference to a peak flag
    const Flag& getPeakFlag(const IndexPair& index) const
    {
      return flags_[index.first][index.second];
    }

    /// Returns mutable reference to a peak flag
    Flag& getPeakFlag(const IndexPair& index)
    {
      return flags_[index.first][index.second];
    }

    /// Returns the default parameters for the algorithm with name @p algorithm_name
    Param getParameters(const String& algorithm_name) const;

protected:

    /// Container for flags attached to input data
    std::vector<std::vector<Flag> > flags_;

  }; // class FeatureFinder

} // namespace OpenMS

