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

#ifndef OPENMS_SIMULATION_EGHFITTER1D_H
#define OPENMS_SIMULATION_EGHFITTER1D_H

#define DEBUG_FEATUREFINDER

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>

namespace OpenMS
{
  /**
    @brief Exponential-Gaussian hybrid distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization.

    @htmlinclude OpenMS_EGHFitter1D.parameters
  */
  class OPENMS_DLLAPI EGHFitter1D :
    public LevMarqFitter1D
  {
public:
    /// Default constructor
    EGHFitter1D();

    /// copy constructor
    EGHFitter1D(const EGHFitter1D & source);

    /// destructor
    ~EGHFitter1D() override;

    /// assignment operator
    virtual EGHFitter1D & operator=(const EGHFitter1D & source);

    /// create new EGHFitter1D object (function needed by Factory)
    static Fitter1D * create()
    {
      return new EGHFitter1D();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "EGHFitter1D";
    }

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, InterpolationModel * & model) override;

protected:

    /// Helper struct (contains the size of an area and a raw data container)
    struct Data
    {
      typedef Peak1D PeakType;
      typedef std::vector<PeakType> RawDataArrayType;

      Size n;
      RawDataArrayType set;
    };

    class EGHFitterFunctor : public LevMarqFitter1D::GenericFunctor
    {
      public:
      EGHFitterFunctor(int dimensions, const EGHFitter1D::Data * data)
      : LevMarqFitter1D::GenericFunctor(dimensions, data->n), m_data(data) {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) override;
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J) override;

      protected:
        const EGHFitter1D::Data * m_data;
    };

    /// Compute start parameter
    virtual void setInitialParameters_(const RawDataArrayType & set);

    /// Parameter of egh - peak height
    CoordinateType height_;
    /// Parameter of egh - tau
    CoordinateType tau_;
    /// Parameter of egh - sigma-square
    CoordinateType sigma_square_;
    /// Parameter of egh - peak retention time
    CoordinateType retention_;

    void updateMembers_() override;
  };

}

#endif // OPENMS_SIMULATION_EGHFITTER1D_H
