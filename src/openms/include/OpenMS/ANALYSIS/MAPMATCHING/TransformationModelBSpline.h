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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELBSPLINE_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELBSPLINE_H

#include <OpenMS/config.h> // is this needed?

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>

namespace OpenMS
{

  /**
    @brief B-spline (non-linear) model for transformations

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelBSpline :
    public TransformationModel
  {
public:
    /**
      @brief Constructor

      @exception Exception::IllegalArgument is thrown if a parameter is invalid.
      @exception Exception::UnableToFit is thrown if the B-spline fit fails.
    */
    TransformationModelBSpline(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelBSpline() override;

    /// Evaluates the model at the given value
    double evaluate(double value) const override;

    using TransformationModel::getParameters;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

protected:
    /// Pointer to the actual B-spline
    BSpline2d* spline_;

    /// Min./max. x value (endpoints of the data range)
    double xmin_, xmax_;

    /// Method to use for extrapolation (beyond 'xmin_'/'xmax_')
    enum { EX_LINEAR, EX_BSPLINE, EX_CONSTANT, EX_GLOBAL_LINEAR } extrapolate_;

    /// Parameters for constant or linear extrapolation 
    double offset_min_, offset_max_, slope_min_, slope_max_;
  };
} // namespace

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELBSPLINE_H
