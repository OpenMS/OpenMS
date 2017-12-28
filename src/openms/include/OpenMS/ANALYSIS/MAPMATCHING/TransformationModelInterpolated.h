// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELINTERPOLATED_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELINTERPOLATED_H

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
  };

} // namespace

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELINTERPOLATED_H
