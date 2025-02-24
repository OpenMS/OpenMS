// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <vector>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /**
      @brief Calculates ticks for a given value range.

      It has only static methods, that's why the constructor is private.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI AxisTickCalculator
  {
public:

    /// Typedef for the grid vector
    typedef std::vector<std::vector<double> > GridVector;

    /**
         @brief Returns a GridVector with ticks for linear scales.

         @param x1 minimum value
         @param x2 maximum value
         @param grid the grid_vector to fill
    */
    static void calcGridLines(double x1, double x2, GridVector & grid);

    /**
         @brief Returns a GridVector with ticks for logarithmic scales.

         @param x1 minimum value
         @param x2 maximum value
         @param grid the grid_vector to fill
    */
    static void calcLogGridLines(double x1, double x2, GridVector & grid);

private:

    ///Constructor: only static methods
    AxisTickCalculator();
  };
}
