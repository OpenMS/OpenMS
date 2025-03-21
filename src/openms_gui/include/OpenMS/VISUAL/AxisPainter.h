// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QPaintEvent>
#include <QPainter>


namespace OpenMS
{
  /**
    @brief Draws a coordinate axis.
    It has only static methods, that's why the constructor is private.
    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI AxisPainter
  {
public:
    /// Typedef for the grid vector
    typedef std::vector<std::vector<double> > GridVector;

    /// Where the axis is placed
    enum Alignment
    {
      TOP,
      BOTTOM,
      LEFT,
      RIGHT
    };

    /// Draws an axis
    static void paint(QPainter * painter, QPaintEvent * e, const double & min, const double & max, const GridVector & grid,
                      const Int width, const Int height, const Alignment alignment, const UInt margin,
                      const bool show_legend, const String& legend, const bool shorten_number,
                      const bool is_log, const bool is_inverse_orientation);
private:
    /// Constructor: only static methods
    AxisPainter();

    /// sets @p short_num to a shortened string representation ("123.4 k/M/G") of @p number
    static void getShortenedNumber_(QString& short_num, double number);

    /// Round to 8 significant digits after comma (and apply log scaling)
    static double scale_(double x, bool is_log);
  };
}
