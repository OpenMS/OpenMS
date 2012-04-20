// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_AXISPAINTER_H
#define OPENMS_VISUAL_AXISPAINTER_H

#include <QtGui/QPainter>

#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <QtGui/QPaintEvent>

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
    static void paint(QPainter * painter, QPaintEvent * e, const DoubleReal & min, const DoubleReal & max, const GridVector & grid,
                      const Int width, const Int height, const Alignment alignment, const UInt margin,
                      const bool show_legend, const String legend, const bool shorten_number,
                      const bool is_log, const bool is_inverse_orientation);
private:
    /// Constructor: only static methods
    AxisPainter();

    /// sets @p short_num to a shortened string representation ("123.4 k/M/G") of @p number
    static void getShortenedNumber_(QString & short_num, DoubleReal number);

    /// Scale axis values to correct value (i.e. reverse log, unit conversion)
    static DoubleReal scale_(DoubleReal x, bool is_log);
  };
}
#endif
