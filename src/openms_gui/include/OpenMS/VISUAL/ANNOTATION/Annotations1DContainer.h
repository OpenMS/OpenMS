// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <list>

#include <QtGui/QPen>

class QPoint;
class QObject;
class QRectF;
class QPainter;

namespace OpenMS
{
  class Annotation1DItem;

  /// Container for annotations to content of Plot1DCanvas
  class Annotations1DContainer :
    public std::list<Annotation1DItem *>
  {
public:
    /// Default constructor
    Annotations1DContainer();

    /// Copy constructor
    Annotations1DContainer(const Annotations1DContainer & rhs);

    /// Assignment operator
    Annotations1DContainer & operator=(const Annotations1DContainer & rhs);

    /// Destructor
    virtual ~Annotations1DContainer();

    using Base = std::list<Annotation1DItem *>;

    /// Iterator for the 1D annotations
    using Iterator = Base::iterator;

    /// Const iterator for the 1D annotations
    using ConstIterator = std::list<Annotation1DItem *>::const_iterator;

    /// Type of the Points
    using PointType = DPosition<2>;

    /// Coordinate type
    using CoordinateType = double;

    /** @brief Returns a pointer to the item at @p pos, or 0, if not existent

            If more than one item's bounding box encloses @p pos , the one in the
            foreground is returned.
    */
    Annotation1DItem * getItemAt(const QPoint & pos) const;

    /// Selects the item at @p pos on the canvas, if it exists.
    void selectItemAt(const QPoint & pos) const;

    /// Deselects the item at @p pos on the canvas, if it exists.
    void deselectItemAt(const QPoint & pos) const;

    /// Selects all items
    void selectAll();

    /// Deselects all items
    void deselectAll();

    /// Removes the selected items
    void removeSelectedItems();

    /// Returns the selected items
    std::vector<Annotation1DItem*> getSelectedItems();

    /// Sets the pen_
    void setPen(const QPen & pen);

    /// Returns the pen_
    const QPen & getPen() const;

    /// Sets the selected_pen_
    void setSelectedPen(const QPen & pen);

    /// Returns the selected_pen_
    const QPen & getSelectedPen() const;

  protected:
    /// call delete on all pointers in the container, without modifying the container
    void deleteAllItems_() const;

    /// The pen used to draw items
    QPen pen_;

    /// The pen used to draw selected items
    QPen selected_pen_;
  };

} // namespace

