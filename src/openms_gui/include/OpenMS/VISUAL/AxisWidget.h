// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// QT
#include <QtWidgets>
class QPaintEvent;

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/VISUAL/AxisPainter.h>

namespace OpenMS
{
  /**
      @brief Widget that represents an axis of a graph.

      Additional to ticks and tick values a label e.g. the unit can be displayed.
      It supports both linear and logarithmic scale.

      @image html AxisWidget.png
      The above image shows a horizontal example axis.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI AxisWidget :
    public QWidget
  {
    Q_OBJECT

public:
    ///Type definitions
    //@{
    ///Vector of vector of doubles that defines the grid
    typedef std::vector<std::vector<double> > GridVector;

    /// constructor
    AxisWidget(const AxisPainter::Alignment alignment, const char * legend = "", QWidget * parent = nullptr);

    /// destructor
    ~AxisWidget() override;

    /// sets the margin on the top/right side (default is 0)
    void setMargin(UInt size);

    /// returns the margin
    UInt margin() const;

    /// enable the display of the legend (default true)
    void showLegend(bool show_legend);

    /// returns true if legend is shown
    bool isLegendShown() const;

    /// sets the legend text
    void setLegend(const String & legend);

    /// returns the actual legend text
    const String & getLegend() const;

    /// returns the currently used grid lines
    const GridVector & gridLines() const;

    /// sets the axis to logarithmic scale
    void setLogScale(bool is_log);

    /// returns true if the axis has logarithmic scale
    bool isLogScale() const;

    /// set true to display the axis label in inverse order (left to right or bottom to top)
    void setInverseOrientation(bool inverse_orientation);

    /// returns if the axis label is displayed in inverse order
    bool hasInverseOrientation() const;

    /// set true to allow for shortened numbers (with k/M/G units) on the axis label
    void setAllowShortNumbers(bool short_nums);

    /// returns the minimum value displayed on the axis
    double getAxisMinimum() const;

    /// returns the maximum value displayed on the axis
    double getAxisMaximum() const;

    /// Actual painting takes place here
    void paint(QPainter * painter, QPaintEvent * e);

public slots:

    ///sets min/max of the axis
    void setAxisBounds(double min, double max);

    /// set maximum number of tick levels ('1' or '2', default: '2')
    void setTickLevel(UInt level);

protected:
    /// Vector that defines the position of the ticks/gridlines and the shown values on axis
    GridVector grid_line_;

    /// format of axis scale (linear or logarithmic)
    bool is_log_;

    /// display of legend enabled or not
    bool show_legend_;

    /// Position of the axis (right, left, top, down as defined in ALIGNMENT_ENUM)
    AxisPainter::Alignment alignment_;

    /// true if axis label are displayed in inverse order (left to right or bottom to top)
    bool is_inverse_orientation_;

    /// margin of axis
    UInt margin_;

    /// minimum value on the axis
    double min_;

    /// maximum value on the axis
    double max_;

    /// text/unit on axis
    String legend_;

    /// maximum number of tick levels (default=2)
    UInt tick_level_;

    /// true if k/M/G units can be used
    bool allow_short_numbers_;

    /// Reimplemented Qt event (calls paint with "this")
    void paintEvent(QPaintEvent *) override;
  };
} // namespace OpenMS

