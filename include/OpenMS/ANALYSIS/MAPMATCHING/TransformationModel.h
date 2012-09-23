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
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_interp.h>

namespace OpenMS
{
  /**
       @brief Base class for transformation models

       Implements the identity (no transformation). Parameters and data are ignored.

       @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModel
  {
public:
    /// Coordinate pair
    typedef std::pair<DoubleReal, DoubleReal> DataPoint;
    /// Vector of coordinate pairs
    typedef std::vector<DataPoint> DataPoints;

    /// Constructor
    TransformationModel() {}

    /// Alternative constructor (derived classes should implement this one!)
    TransformationModel(const TransformationModel::DataPoints &,
                        const Param &) :
      params_() {}

    /// Destructor
    virtual ~TransformationModel() {}

    /// Evaluates the model at the given value
    virtual DoubleReal evaluate(const DoubleReal value) const
    {
      return value;
    }

    /// Gets the (actual) parameters
    void getParameters(Param & params) const
    {
      params = params_;
    }

    /// Gets the default parameters
    static void getDefaultParameters(Param & params)
    {
      params.clear();
    }

protected:
    /// Parameters
    Param params_;
  };


  /**
       @brief Linear model for transformations

       The model can be inferred from data or specified using explicit parameters. If data is given, a least squares fit is used to find the model parameters (slope and intercept). Depending on parameter @p symmetric_regression, a normal regression (@e y on @e x) or symmetric regression (@f$ y - x @f$ on @f$ y + x @f$) is performed.

       Without data, the model can be specified by giving the parameters @p slope and @p intercept explicitly.

       @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelLinear :
    public TransformationModel
  {
public:
    /**
        @brief Constructor

        @exception IllegalArgument is thrown if neither data points nor explicit parameters (slope/intercept) are given.
    */
    TransformationModelLinear(const DataPoints & data, const Param & params);

    /// Destructor
    ~TransformationModelLinear();

    /// Evaluates the model at the given value
    virtual DoubleReal evaluate(const DoubleReal value) const;

    using TransformationModel::getParameters;

    /// Gets the "real" parameters
    void getParameters(DoubleReal & slope, DoubleReal & intercept) const;

    /// Gets the default parameters
    static void getDefaultParameters(Param & params);

    /**
         @brief Computes the inverse

         @exception DivisionByZero is thrown if the slope is zero.
    */
    void invert();

protected:
    /// Parameters of the linear model
    DoubleReal slope_, intercept_;
    /// Was the model estimated from data?
    bool data_given_;
    /// Use symmetric regression?
    bool symmetric_;
  };


  /**
       @brief Interpolation model for transformations

       Between the data points, the interpolation uses the neighboring points. Outside the range spanned by the points, we extrapolate using a line through the first and the last point.

       Different types of interpolation (controlled by the parameter @p interpolation_type) are supported: "linear", "polynomial", "cspline", and "akima". Note that the number of required data points may differ between types.

       @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelInterpolated :
    public TransformationModel
  {
public:
    /**
         @brief Constructor

         @exception IllegalArgument is thrown if there are not enough data points or if an unknown interpolation type is given.
    */
    TransformationModelInterpolated(const DataPoints & data,
                                    const Param & params);

    /// Destructor
    ~TransformationModelInterpolated();

    /// Evaluates the model at the given value
    DoubleReal evaluate(const DoubleReal value) const;

    /// Gets the default parameters
    static void getDefaultParameters(Param & params);

protected:
    /// Data coordinates
    std::vector<double> x_, y_;
    /// Number of data points
    size_t size_;
    /// Look-up accelerator
    gsl_interp_accel * acc_;
    /// Interpolation function
    gsl_interp * interp_;
    /// Linear model for extrapolation
    TransformationModelLinear * lm_;
  };


  /**
       @brief B-spline model for transformations

       In the range of the data points, the transformation is evaluated from a cubic smoothing spline fit to the points. The number of breakpoints is given as a parameter (@p num_breakpoints). Outside of this range, linear extrapolation through the last point with the slope of the spline at that point is used.

       Positioning of the breakpoints is controlled by the parameter @p break_positions. Valid choices are "uniform" (equidistant spacing on the data range) and "quantiles" (equal numbers of data points in every interval).

       @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelBSpline :
    public TransformationModel
  {
public:
    /**
         @brief Constructor

         @exception IllegalArgument is thrown if not enough data points are given or if the required parameter @p num_breakpoints is missing.
    */
    TransformationModelBSpline(const DataPoints & data, const Param & params);

    /// Destructor
    ~TransformationModelBSpline();

    /// Evaluates the model at the given value
    DoubleReal evaluate(const DoubleReal value) const;

    /// Gets the default parameters
    static void getDefaultParameters(Param & params);

protected:
    /**
        @brief Finds quantile values

        @param x Data vector to find quantiles in
        @param quantiles Quantiles to find (values between 0 and 1)
        @param results Resulting quantiles (vector must be already allocated to the correct size!)
    */
    void getQuantiles_(const gsl_vector * x, const std::vector<double> &
                       quantiles, gsl_vector * results);

    /// Computes the B-spline fit
    void computeFit_();

    /// Computes the linear extrapolation
    void computeLinear_(const double pos, double & slope, double & offset,
                        double & sd_err);

    /// Vectors for B-spline computation
    gsl_vector * x_, * y_, * w_, * bsplines_, * coeffs_;
    /// Covariance matrix
    gsl_matrix * cov_;
    /// B-spline workspace
    gsl_bspline_workspace * workspace_;
    /// Number of data points and coefficients
    size_t size_, ncoeffs_;
    // First/last breakpoint
    double xmin_, xmax_;
    /// Parameters for linear extrapolation
    double slope_min_, slope_max_, offset_min_, offset_max_;
    /// Fitting errors of linear extrapolation
    double sd_err_left_, sd_err_right_;
  };

} // end of namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
