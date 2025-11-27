// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit, Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/MISC/CommonDefs.h> // for PointXYType

#include <QtCore/QRectF>
#include <QtCore/QString>

class QPainter;

namespace OpenMS
{
  template<int D> class DimMapper;
  class Gravitator;
  class Plot1DCanvas;

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
    /// Destructor
    virtual ~Annotation1DItem();

    /// Returns the current bounding box of this item on the canvas where it has last been drawn
    const QRectF& boundingBox() const;

    /// Returns true if this item is currently selected on the canvas, else false
    bool isSelected() const;

    /// Sets whether this item is currently selected on the canvas or not
    void setSelected(bool selected);

    /// Sets the text of the item
    void setText(const QString & text);

    /// Returns the text of the item
    const QString& getText() const;

    /// open a GUI input field and let the user edit the text
    /// If the text was changed, true is returned; otherwise false.
    bool editText();

    /// Ensures that the item has coordinates within the visible area of the canvas
    virtual void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) = 0;

    /// Draws the item on @p painter
    virtual void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) = 0;

    /// Moves the item on the drawing canvas; behavior depends on item type and is implemented in the subclasses
    virtual void move(const PointXYType delta, const Gravitator& gr, const DimMapper<2>& dim_mapper) = 0;

    /// Creates a copy of the item on the heap and returns a pointer
    virtual Annotation1DItem* clone() const = 0;

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
