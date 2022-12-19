// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
