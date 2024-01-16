// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
  /**
      @brief A map feature grouping algorithm for labeling techniques with two labels.

      It takes one maps and searches for corresponding features with a defined distance in RT and m/z.

    @htmlinclude OpenMS_FeatureGroupingAlgorithmLabeled.parameters

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmLabeled :
    public FeatureGroupingAlgorithm
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithmLabeled();

    /// Destructor
    ~FeatureGroupingAlgorithmLabeled() override;

    /**
        @brief Applies the algorithm

        @note Exactly one @em input map has to be provided.
        @note The @em output map has to have two file descriptions, containing
        the same file name. The file descriptions have to be labeled 'heavy' and 'light'.

        @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void group(const std::vector<FeatureMap > & maps, ConsensusMap & out) override;

    ///Creates a new instance of this class (for Factory)
    static FeatureGroupingAlgorithm * create()
    {
      return new FeatureGroupingAlgorithmLabeled();
    }

    ///Returns the product name (for the Factory)
    static String getProductName()
    {
      return "labeled";
    }

private:

    ///Copy constructor is not implemented -> private
    FeatureGroupingAlgorithmLabeled(const FeatureGroupingAlgorithmLabeled &);
    ///Assignment operator is not implemented -> private
    FeatureGroupingAlgorithmLabeled & operator=(const FeatureGroupingAlgorithmLabeled &);

  };

} // namespace OpenMS

