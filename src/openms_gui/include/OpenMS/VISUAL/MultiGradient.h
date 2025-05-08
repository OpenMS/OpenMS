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

    /// Returns the default gradient for logarithmic intensity mode
    static MultiGradient getDefaultGradientLogarithmicIntensityMode();

    /// Interpolation mode.
    enum InterpolationMode
    {
      IM_LINEAR,        ///< IM_LINEAR returns the linear interpolation (default).
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
    void insert(double position, QColor color);
    /// removes the color at position @p position
    bool remove(double position);
    /// returns if a value for position @p position exists
    bool exists(double position);
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
    QColor interpolatedColorAt(double position) const;
    /**
        @brief returns the color as @p position with the gradient stretched between @p min and @p max.

        If the @p position is higher or lower than the range [min,max] the highest,
        respectively the lowest, color is returned.
    */
    QColor interpolatedColorAt(double position, double min, double max) const;

    /// activates the precalculation of values (only approximate results are given)
    void activatePrecalculationMode(double min, double max, UInt steps);
    /// deactivates the precalculation of values ( and deletes the precalculated values)
    void deactivatePrecalculationMode();

    /// index of color in precalculated table by position in gradient
    inline Int precalculatedColorIndex( double position ) const
    {
      OPENMS_PRECONDITION(pre_.size() != 0, "MultiGradient::precalculatedColorIndex(double): Precalculation mode not activated!");
      OPENMS_PRECONDITION(position >= pre_min_, (String("MultiGradient::precalculatedColorIndex(double): position ") + position + " out of specified range (" + pre_min_ + "-" + (pre_min_ + pre_size_) + ")!").c_str());

      Int index = (Int)((position - pre_min_) / pre_size_ * pre_steps_);

      return qBound( 0, index, (Int)pre_.size() - 1 );
    }

    /// precalculated color by its index in the table
    inline QColor precalculatedColorByIndex( Int index ) const
    {
      OPENMS_PRECONDITION(pre_.size() != 0, "MultiGradient::precalculatedColorByIndex(Int): Precalculation mode not activated!");
      OPENMS_PRECONDITION( index >= 0, "MultiGradient::precalculatedColorByIndex(Int): negative indexes not allowed");
      OPENMS_PRECONDITION( index < (Int)pre_.size(), (String("MultiGradient::indexedColor(Int): index ") + index + " out of specified range (0-" + pre_.size() + ")!").c_str());

      return pre_[index];
    }

    /**
        @brief Returns a precalculated color.

        If the @p position is out of the range specified in activatePrecalculationMode() the behaviour depends on the debug mode:
        - With debug information an Precondition exception is thrown
        - Without debug information array boundaries are violated, which probably causes a segmentation fault.
    */
    inline QColor precalculatedColorAt(double position) const
    {
      return precalculatedColorByIndex( precalculatedColorIndex( position ) );
    }

    ///return the number of color points
    Size size() const;

    /// size of precalculated colors table
    Size precalculatedSize() const
    {
        return pre_.size();
    }

    /// sets the interpolation mode (default or stairs). Default is linear
    void setInterpolationMode(InterpolationMode mode);
    /// returns the interpolation mode
    InterpolationMode getInterpolationMode() const;

    ///convert to string representation
    std::string toString() const;
    /**
        @brief Sets the gradient by string representation.

        The string representation of a gradient starts with the interpolation mode: "Linear" or "Stairs" and the separator "|".
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
    std::map<double, QColor> pos_col_;
    /// Current interpolation mode
    InterpolationMode interpolation_mode_;
    /// Precalculated colors
    std::vector<QColor> pre_;
    /// Minimum of the precalculated color range
    double pre_min_;
    /// Width of the precalculated color range
    double pre_size_;
    /// Steps of the precalculated color range
    UInt pre_steps_;

  };

}
