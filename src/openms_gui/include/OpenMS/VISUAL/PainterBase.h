// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>
#include <OpenMS/CONCEPT/Types.h>

#include <string_view>

#include <QPainterPath>
#include <QRgb>

class QColor;
class QPainter;
class QPenStyle;
class QPoint;

namespace OpenMS
{
  class String;
  enum class ShapeIcon
  {
    DIAMOND,
    SQUARE,
    CIRCLE,
    TRIANGLE
  };

  /**
   * @brief An empty base class with some static convenience functions
  */
  class OPENMS_GUI_DLLAPI PainterBase
  {
  public:
    /// translates 'diamond', 'square', 'circle', 'triangle' into a ShapeIcon
    /// @throws Exception::InvalidValue otherwise
    static ShapeIcon toShapeIcon(const String& icon);
    
    /// static method to draw a dashed line
    static void drawDashedLine(const QPoint& from, const QPoint& to, QPainter* painter, const QColor& color);

    /// draw a cross at @p position, using a certain size (= width = height) of the cross
    static void drawCross(const QPoint& position, QPainter* painter, const int size = 8);

    /// draw a caret '^' at @p position, using a certain size (= width) of the caret
    static void drawCaret(const QPoint& position, QPainter* painter, const int size = 8);

    /// draw an unfilled diamond at @p position, using a certain size (= width = height) of the diamond
    static void drawDiamond(const QPoint& position, QPainter* painter, const int size = 8);

    /// draws squares, circles etc
    static void drawIcon(const QPoint& pos, const QRgb& color, const ShapeIcon icon, Size s, QPainter& p);

    /// An arrow head which is open, i.e. '>'
    static QPainterPath getOpenArrow(int arrow_width);
    /// An arrow head which is closed, i.e. a triangle
    static QPainterPath getClosedArrow(int arrow_width);

    /**
     * \brief 
     * \param painter The painter to paint with
     * \param pen For setting line width and color
     * \param start Start position of the line
     * \param end End position of the line
     * \param arrow_start An (optional) arrow head. Use 'getOpenArrow' or 'getClosedArrow' for predefined arrows
     * \param arrow_end  An (optional) arrow tail. Use 'getOpenArrow' or 'getClosedArrow' for predefined arrows
     * \return The bounding rectangle of the line and arrows (if any)
     */
    static QRectF drawLineWithArrows(QPainter* painter, const QPen& pen, const QPoint& start, const QPoint& end, 
                                     const QPainterPath& arrow_start = QPainterPath(),
                                     const QPainterPath& arrow_end = QPainterPath());
  };

} // namespace OpenMS
