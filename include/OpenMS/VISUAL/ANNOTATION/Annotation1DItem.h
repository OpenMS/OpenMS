// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ANNOTATION_ANNOTATION1DITEM_H
#define OPENMS_VISUAL_ANNOTATION_ANNOTATION1DITEM_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <QtCore/QRectF>
#include <QtCore/QString>

class QPainter;

namespace OpenMS
{
  class Spectrum1DCanvas;

  /** @brief An abstract class acting as an interface for the different 1D annotation items.

            This is an abstract polymorphic class which acts as an interface between its
            subclasses and all containers and methods that contain or handle Annotation1DItem
            objects.

            If you want to add new kinds of annotation items, inherit this class,
            implement the pure virtual methods, and add everything else the item should
            have or be capable of.

    */
  class Annotation1DItem
  {

public:
    /// Type of the Points
    typedef DPosition<2> PointType;

    /// Intensity type
    typedef Real IntensityType;

    /// Coordinate type
    typedef DoubleReal CoordinateType;

    /// Destructor
    virtual ~Annotation1DItem();

    /// Returns the current bounding box of this item on the canvas where it has last been drawn
    const QRectF & boundingBox() const;

    /// Returns true if this item is currently selected on the canvas, else false
    bool isSelected() const;

    /// Sets whether this item is currently selected on the canvas or not
    void setSelected(bool selected);

    /// Sets the text of the item
    void setText(const QString & text);

    /// Returns the text of the item
    const QString & getText() const;

    /// Ensures that the item has coordinates within the visible area of the canvas
    virtual void ensureWithinDataRange(Spectrum1DCanvas * const canvas) = 0;

    /// Draws the item on @p painter
    virtual void draw(Spectrum1DCanvas * const canvas, QPainter & painter, bool flipped = false) = 0;

    /// Moves the item; behaviour depends on item type and is implemented in the subclasses
    virtual void move(const PointType & delta) = 0;

protected:
    /// Constructor
    Annotation1DItem(const QString & text);

    /// Copy constructor
    Annotation1DItem(const Annotation1DItem & rhs);

    /// Draws the bounding_box_
    void drawBoundingBox_(QPainter & painter);

    /// The current bounding box of this item on the canvas where it has last been drawn
    QRectF bounding_box_;

    /// Determines whether this item is currently selected on the canvas
    bool selected_;

    /// The displayed text
    QString text_;
  };
} // namespace OpenMS

#endif
