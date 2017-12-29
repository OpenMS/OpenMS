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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_EGHMODEL_H
#define OPENMS_SIMULATION_EGHMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <boost/math/tr1.hpp>

namespace OpenMS
{
  /**
      @brief Exponential-Gaussian hybrid distribution model for elution profiles.

  Lan K, Jorgenson JW.
  A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks.
  Journal of Chromatography A. 2001;915(1-2):1-13.
  Available at: http://linkinghub.elsevier.com/retrieve/pii/S0021967301005945

  @htmlinclude OpenMS_EGHModel.parameters

  */
  class OPENMS_DLLAPI EGHModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef Math::BasicStatistics<CoordinateType> BasicStatistics;
    typedef LinearInterpolation::container_type ContainerType;

    /// Default constructor
    EGHModel();

    /// copy constructor
    EGHModel(const EGHModel & source);

    /// destructor
    ~EGHModel() override;

    /// assignment operator
    virtual EGHModel & operator=(const EGHModel & source);

    /// create new ElutionModel object (needed by Factory)
    static BaseModel<1> * create()
    {
      return new EGHModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "EGHModel";
    }

    /// set offset without being computing all over and without any discrepancy
    void setOffset(CoordinateType offset) override;

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /// get the center of the Gaussian model i.e. the position of the maximum
    CoordinateType getCenter() const override;

protected:
    CoordinateType  min_;
    CoordinateType  max_;
    BasicStatistics statistics_;
    CoordinateType  height_;     // H in paper
    CoordinateType  apex_rt_;

    CoordinateType  A_;
    CoordinateType  B_;

    CoordinateType  tau_;
    CoordinateType  sigma_square_;
    CoordinateType  sigma_square_2_;


    void updateMembers_() override;

    /// Computes a left & right boundary for the EGH Profile and sets the internal parameters accordingly
    void computeBoundaries_();

    /**
     * @brief Evaluate the EGH function at position rt
     *
     * @param rt        The position where the EGH function should be evaluated. Note that this is the position without the RT offset, meaning that the EGH apex is at position 0
     * @param egh_value The computed value
     */
    inline void evaluateEGH_(CoordinateType & rt, CoordinateType & egh_value)
    {
      CoordinateType denominator = sigma_square_2_ + tau_ * rt;

      if (denominator > 0)
      {
        // evaluate egh ->
        egh_value = height_ * exp(
          (-1 * rt * rt) / denominator
          );
      }
      else
      {
        egh_value = 0.0;
      }
    }

  };

} // namespace OpenMS

#endif // OPENMS_SIMULATION_ELUTIONMODEL_H
