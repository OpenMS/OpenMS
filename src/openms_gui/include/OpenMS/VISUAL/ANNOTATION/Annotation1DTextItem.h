// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Johannes Junker, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

namespace OpenMS
{
  /** @brief An annotation item which represents an arbitrary text on the canvas.
          @see Annotation1DItem
  */
  class Annotation1DTextItem :
    public Annotation1DItem
  {
public:

    /// Constructor
    Annotation1DTextItem(const PointXYType& position, const QString& text, const int flags = Qt::AlignCenter)
      : Annotation1DItem(text), position_(position), flags_(flags)
    {
    }
    /// Copy constructor
    Annotation1DTextItem(const Annotation1DTextItem & rhs) = default;

    /// Destructor
    ~Annotation1DTextItem() override = default;

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override;

    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override;

    // Docu in base class
    void move(const PointXYType delta, const Gravitator& gr, const DimMapper<2>& dim_mapper) override;

    /// Sets the position of the item (in X-Y coordinates)
    void setPosition(const PointXYType& position)
    {
      position_ = position;
    }

    /// Returns the position of the item (in X-Y coordinates)
    const PointXYType& getPosition() const
    {
      return position_;
    }

    /// Set Qt flags (default: Qt::AlignCenter)
    void setFlags(int flags)
    {
      flags_ = flags;
    }

    /// Get Qt flags
    int getFlags() const
    {
      return flags_;
    }

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DTextItem(*this);
    }

  protected:
    /// The position of the item as a datatype, e.g. Peak1D
    PointXYType position_;

    int flags_;
  };
} // namespace OpenMS

