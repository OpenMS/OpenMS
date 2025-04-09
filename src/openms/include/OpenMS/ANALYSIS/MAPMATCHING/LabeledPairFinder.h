// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#include <cmath>

namespace OpenMS
{
  /**
      @brief The LabeledPairFinder allows the matching of labeled features (features with a fixed distance).

      Finds feature pairs that have a defined distance in RT and m/z in the same map.

      @htmlinclude OpenMS_LabeledPairFinder.parameters

      @todo Implement support for labeled MRM experiments, Q1 m/z value and charges. (Andreas)
      @todo Implement support for more than one mass delta, e.g. from missed cleavages and so on (Andreas)

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI LabeledPairFinder :
    public BaseGroupFinder
  {

public:

    /// Default constructor
    LabeledPairFinder();

    /// Destructor
    inline ~LabeledPairFinder() override
    {
    }

    /**
      @brief Run the algorithm

      @note Exactly one @em input map has to be provided.
      @note The @em output map has to have two file descriptions, containing
      the same file name. The file descriptions have to be labeled 'heavy'
      and 'light'.

      @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<ConsensusMap> & input_maps, ConsensusMap & result_map) override;

protected:

    /// return the p-value at position x for the bi-Gaussian distribution with mean @p m and standard deviation @p sig1 (left) and @p sig2 (right)
    inline double PValue_(double x, double m, double sig1, double sig2)
    {
      if (m < x)
      {
        return 1 - std::erf((x - m) / sig2 / 0.707106781);
      }
      else
      {
        return 1 - std::erf((m - x) / sig1 / 0.707106781);
      }
    }

private:

    /// Copy constructor not implemented => private
    LabeledPairFinder(const LabeledPairFinder & source);

    /// Assignment operator not implemented => private
    LabeledPairFinder & operator=(const LabeledPairFinder & source);

  };   // end of class LabeledPairFinder

} // end of namespace OpenMS

