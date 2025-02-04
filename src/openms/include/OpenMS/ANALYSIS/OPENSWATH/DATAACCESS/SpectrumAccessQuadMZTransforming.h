// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>

namespace OpenMS
{
  /**
   * @brief A transforming m/z wrapper around spectrum access using a quadratic equation. 
   *
   * For each spectrum access, each m/z value is transformed using the equation 
   *    mz = a + b * mz + c * mz^2
   * 
   * This can be used to implement an on-line mass correction for TOF
   * instruments (for example).
   *
   */
  class OPENMS_DLLAPI SpectrumAccessQuadMZTransforming :
    public SpectrumAccessTransforming
  {
public:

    /** @brief Constructor
     * @param sptr The spectrum to work on
     * @param a Regression parameter 0
     * @param b Regression parameter 1
     * @param c Regression parameter 2
     * @param ppm Whether the transformation should be applied in ppm domain
     *            (if false, it is applied directly in m/z domain)
     *
    */
    explicit SpectrumAccessQuadMZTransforming(OpenSwath::SpectrumAccessPtr sptr,
        double a, double b, double c, bool ppm);
        
    ~SpectrumAccessQuadMZTransforming() override;

    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

private:

    double a_;
    double b_;
    double c_;
    bool ppm_;

  };
}

