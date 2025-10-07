// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <iosfwd>
#include <map>

namespace OpenMS
{
  /**
    @brief Generic description of a coordinate transformation.

    This description primarily stores data points (coordinate pairs) from which a @ref TransformationModel "transformation model" can be estimated. Applying the transformation to a coordinate (via @p apply) then means evaluating the model at that coordinate.

    The following models are available:
    - @p none (TransformationModel): \f$ f(x) = x \f$ (identity)
    - @p identity: Same as @p none, but intended for reference files (used to indicate that no other model should be fit, because the identity is already optimal).
    - @p linear (TransformationModelLinear): \f$ f(x) = slope * x + intercept \f$
    - @p interpolated (TransformationModelInterpolated): Interpolation between pairs, extrapolation using first and last pair. Supports different interpolation types.
    - @p b-spline (TransformationModelBSpline): Non-linear smoothing spline, with different options for extrapolation.
    - @b lowess (TransformationModelLowess): Non-linear smoothing via local regression, with different options for extrapolation.

    @remark TransformationDescription stores data points, TransformationModel stores parameters. That way, data can be modeled using different models/parameters, and models can still keep a representation of the data in the format they need (if at all).

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationDescription
  {
    // friend class MapAlignmentAlgorithm;

public:

    /** @brief Summary statistics before/after applying the transformation.
              For deviations before/after transformation, the percentiles 
              100, 99, 95, 90, 75, 50, 25 are returned.
     */ 
    struct TransformationStatistics
    {
      // default constructor
      TransformationStatistics() = default;

      // copy constructor
      TransformationStatistics(const TransformationStatistics& rhs) = default;

      // copy assignment
      TransformationStatistics& operator=(const TransformationStatistics& rhs) = default;

      std::vector<Size> percents = {100, 99, 95, 90, 75, 50, 25};  // TODO: use constexpr array
      double xmin = 0; ///< smallest x value before transformation
      double xmax = 0; ///< largest x value before transformation
      double ymin = 0; ///< smallest y value before transformation
      double ymax = 0; ///< largest y value before transformation
      
      std::map<Size, double> percentiles_before; ///< percentiles of x/y deviations before transformation
      std::map<Size, double> percentiles_after; ///< percentiles of x/y deviations after transformation
    };
    

    /// Coordinate pair
    typedef TransformationModel::DataPoint DataPoint;
    /// Vector of coordinate pairs
    typedef TransformationModel::DataPoints DataPoints;

    /// Default constructor
    TransformationDescription();
    
    /// Constructor from data
    explicit TransformationDescription(const DataPoints& data);
    
    /// Destructor
    ~TransformationDescription();

    /// Copy constructor
    TransformationDescription(const TransformationDescription& rhs);
    /// Assignment operator
    TransformationDescription& operator=(const TransformationDescription& rhs);

    /// Fits a model to the data
    void fitModel(const String& model_type, const Param& params = Param());

    /**
      @brief Applies the transformation to @p value.

      Returns the result of evaluating the fitted model at @p value.
      Returns @p value unchanged if no model was fitted.
    */
    double apply(double value) const;

    /// Gets the type of the fitted model
    const String& getModelType() const;

    /// Gets the possible types of models
    static void getModelTypes(StringList& result);

    /**
      @brief Sets the data points

      Removes the model that was previously fitted to the data (if any).
    */
    void setDataPoints(const DataPoints& data);

    /**
      @brief Sets the data points (backwards-compatible overload)

      Removes the model that was previously fitted to the data (if any).
    */
    void setDataPoints(const std::vector<std::pair<double, double> >& data);

    /// Returns the data points
    const DataPoints& getDataPoints() const;

    /// Non-mutable access to the model parameters
    const Param& getModelParameters() const;

    /// Computes an (approximate) inverse of the transformation
    void invert();

    /**
       @brief Get the deviations between the data pairs

       @param diffs Output
       @param do_apply Get deviations after applying the model?
       @param do_sort Sort @p diffs before returning?
    */
    void getDeviations(std::vector<double>& diffs, bool do_apply = false,
                       bool do_sort = true) const;

    /// Get summary statistics (ranges and errors before/after)
    TransformationStatistics getStatistics() const;

    /// Print summary statistics for the transformation
    void printSummary(std::ostream& os) const;

protected:
    /// Data points
    DataPoints data_;
    /// Type of model
    String model_type_;
    /// Pointer to model
    TransformationModel* model_;
  };

} // end of namespace OpenMS

