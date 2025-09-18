// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>

namespace OpenMS
{
  /**
    @brief Interpolation model for transformations

    Between the data points, the interpolation uses the neighboring points to
    interpolate. The following interpolation methods are available:

    - linear: Linearly interpolate between neighboring points
    - cspline: Use a cubic spline to interpolate between neighboring points
    - akima: Use an akima spline to interpolate between neighboring points (less affected by outliers)

    Outside the range spanned by the points, we extrapolate using one of the following methods:

    - two-point-linear: Uses a line through the first and last point to extrapolate 
    - four-point-linear: Uses a line through the first and second point to
                         extrapolate in front and and a line through the last
                         and second-to-last point in the end. If the data is
                         non-linear, this may yield better approximations for
                         extrapolation.
    - global-linear: Uses a linear regression to fit a line through all data
                     points and use it for extrapolation. Note that
                     global-linear extrapolation may \b not be continuous with
                     the interpolation at the border.

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelInterpolated :
    public TransformationModel
  {
public:

    /**
      @brief Constructor

      @param data The known data points.
      @param params Param object holding information on which model to choose.

      @exception IllegalArgument is thrown if there are not enough data points or if an unknown interpolation type is given.
    */
    TransformationModelInterpolated(const DataPoints& data, const Param& params);
    TransformationModelInterpolated(const std::vector<std::pair<double,double>>& data, const Param& params, bool preprocess);

    /// Destructor
    ~TransformationModelInterpolated() override;

    /**
     * @brief Evaluate the interpolation model at the given value
     *
     * @param value The position where the interpolation should be evaluated.
     *
     * @return The interpolated value.
     */
    double evaluate(double value) const override;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

    /**
     * @brief The class defines a generic interpolation technique used in the
     * TransformationModelInterpolated.
     *
     * @note This class is nested in cpp file as we don't want this to be part
     *       of the public interface nor to be exposed to derived or other classes.
     */
    class Interpolator
    {
public:

      /**
       * @brief Initialize the Interpolator.
       *
       * @param x The x data.
       * @param y The y data.
       */
      virtual void init(std::vector<double>& x, std::vector<double>& y) = 0;

      /**
       * @brief Evaluate the underlying interpolation at a specific position x.
       *
       * @param x The position where the interpolation should be evaluated.
       *
       * @return The interpolated value.
       */
      virtual double eval(const double& x) const = 0;

      /**
       * @brief d'tor.
       */
      virtual ~Interpolator() {}
    };

private:
    /// Data coordinates x
    std::vector<double> x_;

    /// Data coordinates y
    std::vector<double> y_;

    /// Interpolation function
    Interpolator* interp_;

    /// Linear model for extrapolation (front)
    TransformationModelLinear* lm_front_;

    /// Linear model for extrapolation (back)
    TransformationModelLinear* lm_back_;

    /// Preprocesses the incoming data and fills the (private) vectors x_ and y_
    void preprocessDataPoints_(const DataPoints& data);

    /// Preprocesses the incoming data and fills the (private) vectors x_ and y_
    void preprocessDataPoints_(const std::vector<std::pair<double,double>>& data);
  };

} // namespace

