// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>

namespace OpenMS
{
  /**
  @brief Generic description of a coordinate transformation.

  This description primarily stores data points (coordinate pairs) from which a @ref TransformationModel "transformation model" can be estimated. Applying the transformation to a coordinate (via @p apply) then means evaluating the model at that coordinate.

  The following models are available:
  - @p none (TransformationModel): \f$ f(x) = x \f$ (identity)
  - @p identity: Same as @p none, but intended for reference files (used to indicate that no other model should be fit, because the identity is already optimal).
  - @p linear (TransformationModelLinear): \f$ f(x) = slope * x + intercept \f$
  - @p interpolated (TransformationModelInterpolated): Smoothing cubic B-spline.

  @remark TransformationDescription stores data points, TransformationModel stores parameters. That way, data can be modeled using different models/parameters, and models can still keep a representation of the data in the format they need (if at all).

  @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationDescription
  {
    // friend class MapAlignmentAlgorithm;

public:

    /// Coordinate pair
    typedef TransformationModel::DataPoint DataPoint;
    /// Vector of coordinate pairs
    typedef TransformationModel::DataPoints DataPoints;

    /// Default constructor
    TransformationDescription();
    /// Constructor from data
    explicit TransformationDescription(const DataPoints & data);
    /// Destructor
    ~TransformationDescription();

    /// Copy constructor
    TransformationDescription(const TransformationDescription & rhs);
    /// Assignment operator
    TransformationDescription & operator=(const TransformationDescription & rhs);

    /// Fits a model to the data
    void fitModel(const String & model_type, const Param & params = Param());

    /**
         @brief Applies the transformation to @p value.

         Returns the result of evaluating the fitted model at @p value.
         Returns @p value unchanged if no model was fitted.
    */
    DoubleReal apply(DoubleReal value) const;

    /// Gets the type of the fitted model
    const String & getModelType() const;

    /// Gets the possible types of models
    static void getModelTypes(StringList & result);

    /**
         @brief Sets the data points

         Removes the model that was previously fitted to the data (if any).
    */
    void setDataPoints(const DataPoints & data);

    /// Returns the data points
    const DataPoints & getDataPoints() const;

    /// Non-mutable access to the model parameters
    const Param& getModelParameters() const;

    /// Computes an (approximate) inverse of the transformation
    void invert();

protected:
    /// Data points
    DataPoints data_;
    /// Type of model
    String model_type_;
    /// Pointer to model
    TransformationModel * model_;
  };

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H
