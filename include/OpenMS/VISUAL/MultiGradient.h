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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MULTIGRADIENT_H
#define OPENMS_VISUAL_MULTIGRADIENT_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

//QT
#include <QtGui/QColor>

//STL
#include <map>
#include <vector>
#include <cmath>

namespace OpenMS
{

  /**
      @brief A gradient of multiple colors and arbitrary distances between colors.

      Positions associated with numbers range from 0 to 100.
      There is always a color associated with position 0 and 100.
      Stretching the gradient to a specified range, and precalculation and
      caching is also possible.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI MultiGradient
  {
public:
    /// Returns the default gradient for linear intensity mode
    static MultiGradient getDefaultGradientLinearIntensityMode();

    /// Returns the default gradient for linear intensity mode
    static MultiGradient getDefaultGradientLogarithmicIntensityMode();

    /// Interploation mode.
    enum InterpolationMode
    {
      IM_LINEAR,        ///< IM_LINEAR returns the linear interploation (default).
      IM_STAIRS         ///< IM_STAIRS returns the color of the next lower position
    };

    /// Constructor
    MultiGradient();

    /// Copy constructor
    MultiGradient(const MultiGradient & multigradient);

    /// Destructor
    ~MultiGradient();

    /// Assignment operator
    MultiGradient & operator=(const MultiGradient & rhs);

    /// sets or replaces the color at position @p position
    void insert(DoubleReal position, QColor color);
    /// removes the color at position @p position
    bool remove(DoubleReal position);
    /// returns if a value for position @p position exists
    bool exists(DoubleReal position);
    /**
      @brief returns the position of the @p index -th point

      @exception Exception::IndexOverflow is thrown for a too large index
    */
    UInt position(UInt index);
    /**
      @brief returns the color of the @p index -th point

      @exception Exception::IndexOverflow is thrown for a too large index
    */
    QColor color(UInt index);


    /**
        @brief Returns the color as @p position.

        If the @p position is higher or lower than the range [0,100] the highest,
        respectively the lowest, color is returned.
    */
    QColor interpolatedColorAt(DoubleReal position) const;
    /**
        @brief returns the color as @p position with the gradient stretched between @p min and @p max.

        If the @p position is higher or lower than the range [min,max] the highest,
        respectively the lowest, color is returned.
    */
    QColor interpolatedColorAt(DoubleReal position, DoubleReal min, DoubleReal max) const;

    /// activates the precalculation of values (only approximate results are given)
    void activatePrecalculationMode(DoubleReal min, DoubleReal max, UInt steps);
    /// deactivates the precalculation of values ( and deletes the precalculated values)
    void deactivatePrecalculationMode();
    /**
        @brief Returns a precalculated color.

        If the @p position is out of the the range specified in activatePrecalculationMode() the behaviour depends on the debug mode:
        - With debug information an Precondition exception is thrown
        - Wihtout debug information array boundaries are violated, which probably causes a segmentation fault.
    */
    inline QColor precalculatedColorAt(DoubleReal position) const
    {
      OPENMS_PRECONDITION(pre_.size() != 0, "MultiGradient::precalculatedColorAt(DoubleReal): Precalculation mode not activated!");
      OPENMS_PRECONDITION(position >= pre_min_, (String("MultiGradient::precalculatedColorAt(DoubleReal): position ") + position + " out of specified range (" + pre_min_ + "-" + (pre_min_ + pre_size_) + ")!").c_str());
      OPENMS_PRECONDITION(position <= pre_min_ + pre_size_ + std::numeric_limits<DoubleReal>::epsilon() * (pre_min_ + pre_size_), (String("MultiGradient::precalculatedColorAt(DoubleReal): position ") + position + " out of specified range (" + pre_min_ + "-" + (pre_min_ + pre_size_) + ")!").c_str());

      return pre_[(UInt)((position - pre_min_) / pre_size_ * pre_steps_)];
    }

    ///return the number of color points
    Size size() const;

    /// sets the interploation mode (default or stairs). Default is linear
    void setInterpolationMode(InterpolationMode mode);
    /// returns the interpolaion mode
    InterpolationMode getInterpolationMode() const;

    ///convert to string representation
    std::string toString() const;
    /**
        @brief Sets the gradient by string representation.

        The string represenation of a gradient starts with the interpolation mode: "Linear" or "Stairs" and the separator "|".
        It is followed by an arbitrary number of integer-color-pairs.

  Such a pair consists of floating point number (0.0-100.0) followed by a comma and
        a "#". Then follows a color in RGB notation "#RRGGBB" and finally a semicolon.

        Examples are:
        <UL>
            <LI> "Linear|0,#ffff00;100,#000000"
    <LI> "Stairs|0,#ffff00;11.5,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000"
        </UL>
    */
    void fromString(const std::string & gradient);

protected:
    /// Map of index and color
    std::map<DoubleReal, QColor> pos_col_;
    /// Current interpolation mode
    InterpolationMode interpolation_mode_;
    /// Precalculated colors
    std::vector<QColor> pre_;
    /// Minimum of the precalculated color range
    DoubleReal pre_min_;
    /// Width of the precalculated color range
    DoubleReal pre_size_;
    /// Steps of the precalculated color range
    UInt pre_steps_;

  };

}
#endif // OPENMS_VISUAL_MULTIGRADIENT_H
