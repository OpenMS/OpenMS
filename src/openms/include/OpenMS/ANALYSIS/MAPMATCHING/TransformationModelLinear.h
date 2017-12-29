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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELLINEAR_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELLINEAR_H

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

  /**
    @brief Linear model for transformations

    The model can be inferred from data or specified using explicit parameters. 
    If data is given, a least squares fit is used to find the model parameters (slope and intercept). 
    Depending on parameter @p symmetric_regression, a normal regression (@e y on @e x) or
    symmetric regression (@f$ y - x @f$ on @f$ y + x @f$) is performed.

    Without data, the model can be specified by giving the parameters @p slope, @p intercept, 
    @p x_weight, @p y_weight explicitly.

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
    TransformationModelLinear(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelLinear() override;

    /// Evaluates the model at the given value
    double evaluate(double value) const override;

    using TransformationModel::getParameters;

    /// Gets the "real" parameters
    void getParameters(double& slope, double& intercept, String& x_weight, String& y_weight, double& x_datum_min, double& x_datum_max, double& y_datum_min, double& y_datum_max) const;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

    /**
     @brief Computes the inverse

     @exception DivisionByZero is thrown if the slope is zero.
    */
    void invert();

protected:
    /// Parameters of the linear model
    double slope_, intercept_;
    /// Was the model estimated from data?
    bool data_given_;
    /// Use symmetric regression?
    bool symmetric_;
  };
} // namespace

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODELLINEAR_H
