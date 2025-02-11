// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

namespace OpenMS
{
  class EmpiricalFormula;

  /**
        @brief Isotope distribution approximated using linear interpolation.

    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query).
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine

    Peak widening is achieved by either a Gaussian or Lorentzian shape.

        @htmlinclude OpenMS_IsotopeModel.parameters
    */
  class OPENMS_DLLAPI IsotopeModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef InterpolationModel::CoordinateType IntensityType;

    enum Averagines {C = 0, H, N, O, S, AVERAGINE_NUM};

    /// Default constructor
    IsotopeModel();

    ///  copy constructor
    IsotopeModel(const IsotopeModel & source);

    /// destructor
    ~IsotopeModel() override;

    /// assignment operator
    virtual IsotopeModel & operator=(const IsotopeModel & source);

    UInt getCharge() const;

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over.
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model.
    */
    void setOffset(CoordinateType offset) override;

    CoordinateType getOffset();

    /// return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
    EmpiricalFormula getFormula();

    /// set sample/supporting points of interpolation
    virtual void setSamples(const EmpiricalFormula & formula);

    /// set sample/supporting points of interpolation (from base class)
    using InterpolationModel::setSamples;

    /** @brief get the center of the Isotope model

         This is a m/z-value not necessarily the monoisotopic mass.
    */
    CoordinateType getCenter() const override;

    /** @brief the Isotope distribution (without widening) from the last setSamples() call

      Useful to determine the number of isotopes that the model contains and their position

    */
    const IsotopeDistribution & getIsotopeDistribution() const;


protected:
    CoordinateType isotope_stdev_;
    CoordinateType isotope_lorentz_fwhm_;

    UInt charge_;
    CoordinateType mean_;
    CoordinateType monoisotopic_mz_;
    double averagine_[AVERAGINE_NUM];
    Int max_isotope_;
    double trim_right_cutoff_;
    double isotope_distance_;
    IsotopeDistribution isotope_distribution_;

    void updateMembers_() override;

  };
}

