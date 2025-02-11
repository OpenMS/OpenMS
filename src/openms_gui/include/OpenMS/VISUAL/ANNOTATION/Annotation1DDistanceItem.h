// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit, Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <vector>

namespace OpenMS
{
  /** @brief An annotation item which represents a measured distance between two peaks.
      @see Annotation1DItem
  */
  class Annotation1DDistanceItem :
    public Annotation1DItem
  {

public:
    /**
     * \brief 
     * \param text The text to display between the two points
     * \param start_point Start point in XY unit coordinates
     * \param end_point End point in XY unit coordinates 
     * \param swap_ends_if_negative Make sure the distance is positive when creating the distance item?
     */
    Annotation1DDistanceItem(const QString& text, const PointXYType& start_point, const PointXYType& end_point, const bool swap_ends_if_negative = true);
    /// Copy constructor
    Annotation1DDistanceItem(const Annotation1DDistanceItem & rhs) = default;
    /// Destructor
    ~Annotation1DDistanceItem() override = default;

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override;

    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override;
    
    // Docu in base class
    void move(const PointXYType delta, const Gravitator& gr, const DimMapper<2>& dim_mapper) override;

    /// Returns the start point
    const PointXYType& getStartPoint() const
    {
      return start_point_;
    }

    /// Returns the end point
    const PointXYType& getEndPoint() const
    {
      return end_point_;
    }

    /**
     * \brief Compute the (negative) euclidean distance between start and endpoint.
     *
     * If the startpoint is closer to (0,0) than the endpoint, the distance will be positive; otherwise negative.
     *
     * \return sign * sqrt(deltaX^2 + deltaY^2), where deltaX/Y is the difference between start and endpoint in dimension X/Y
     */
    double getDistance() const;

    /// Set tick lines for the distance item in unit XY coordinates (the gravity dimension is ignored)
    void setTicks(const std::vector<PointXYType>& ticks);

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DDistanceItem(*this);
    }

  protected:
    /// The start point of the measured distance line (in XY data coordinates)
    PointXYType start_point_;
    /// The end point of the measured distance line (in XY data coordinates)
    PointXYType end_point_;
    /// Additional tick lines for the distance item (the gravity dimension is ignored)
    std::vector<PointXYType> ticks_;

  };
} // namespace OpenMS

