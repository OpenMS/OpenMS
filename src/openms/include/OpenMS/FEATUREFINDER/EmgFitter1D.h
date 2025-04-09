// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
  /**
    @brief Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization.

    @htmlinclude OpenMS_EmgFitter1D.parameters
  */
  class OPENMS_DLLAPI EmgFitter1D :
    public LevMarqFitter1D
  {
public:

    /// Default constructor
    EmgFitter1D();

    /// copy constructor
    EmgFitter1D(const EmgFitter1D& source);

    /// destructor
    ~EmgFitter1D() override;

    /// assignment operator
    virtual EmgFitter1D& operator=(const EmgFitter1D& source);

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType& range, std::unique_ptr<InterpolationModel>& model) override;

protected:
    /// Helper struct (contains the size of an area and a raw data container)
    struct Data
    {
      typedef Peak1D PeakType;
      typedef std::vector<PeakType> RawDataArrayType;

      Size n;
      RawDataArrayType set;
    };

    class EgmFitterFunctor :
      public LevMarqFitter1D::GenericFunctor
    {
public:
      EgmFitterFunctor(int dimensions, const EmgFitter1D::Data* data) :
        LevMarqFitter1D::GenericFunctor(dimensions,
                                        static_cast<int>(data->n)),
        m_data(data)
      {}

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const override;
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const override;

protected:
      const EmgFitter1D::Data* m_data;
      static const EmgFitter1D::CoordinateType c;
      static const EmgFitter1D::CoordinateType sqrt2pi;
      static const EmgFitter1D::CoordinateType emg_const;
      static const EmgFitter1D::CoordinateType sqrt_2;
    };

    /// Compute start parameter
    virtual void setInitialParameters_(const RawDataArrayType& set);
    /// Compute start parameters using method of moments (usually reduces nr. of iterations needed at some
    /// additional one-time costs
    void setInitialParametersMOM_(const RawDataArrayType& set);

    /// Parameter of emg - peak height
    CoordinateType height_;
    /// Parameter of emg - peak width
    CoordinateType width_;
    /// Parameter of emg - peak symmetry
    CoordinateType symmetry_;
    /// Parameter of emg - peak retention time
    CoordinateType retention_;

    void updateMembers_() override;
  };

}

