// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{
  /**
        @brief Extended isotope distribution approximated using linear interpolation.

    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query).
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine

    Peak widening is achieved by a Gaussian shape.

        @htmlinclude OpenMS_ExtendedIsotopeModel.parameters
    */
  class OPENMS_DLLAPI ExtendedIsotopeModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef InterpolationModel::CoordinateType IntensityType;

    enum Averagines {C = 0, H, N, O, S, AVERAGINE_NUM};

    /// Default constructor
    ExtendedIsotopeModel();

    ///  copy constructor
    ExtendedIsotopeModel(const ExtendedIsotopeModel & source);

    /// destructor
    ~ExtendedIsotopeModel() override;

    /// assignment operator
    virtual ExtendedIsotopeModel & operator=(const ExtendedIsotopeModel & source);

    UInt getCharge() const;

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over.
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model.
    */
    void setOffset(CoordinateType offset) override;

    CoordinateType getOffset();

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /** @brief get the monoisotopic mass of the Isotope model
    */
    CoordinateType getCenter() const override;

protected:
    CoordinateType isotope_stdev_;
    UInt charge_;
    CoordinateType monoisotopic_mz_;
    double averagine_[AVERAGINE_NUM];
    Int max_isotope_;
    double trim_right_cutoff_;
    double isotope_distance_;

    void updateMembers_() override;
  };
}

