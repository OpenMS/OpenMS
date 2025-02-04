// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <memory>

namespace OpenMS
{
  class InterpolationModel;

  /**
  @brief Abstract base class for all 1D-dimensional model fitter.

  @htmlinclude OpenMS_Fitter1D.parameters

  @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI Fitter1D : public DefaultParamHandler
  {
  public:
    /// IndexSet
    typedef IsotopeCluster::IndexSet IndexSet;
    /// IndexSet with charge information
    typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;
    /// Single coordinate
    typedef Feature::CoordinateType CoordinateType;
    /// Quality of a feature
    typedef Feature::QualityType QualityType;
    /// Peak type data point type
    typedef Peak1D PeakType;
    /// Peak type data container type using for the temporary storage of the input data
    typedef std::vector<PeakType> RawDataArrayType;
    /// Peak type data iterator
    typedef RawDataArrayType::iterator PeakIterator;

    /// Default constructor.
    Fitter1D();

    /// copy constructor
    Fitter1D(const Fitter1D& source);

    /// destructor
    ~Fitter1D() override;

    /// assignment operator
    Fitter1D& operator=(const Fitter1D& source);

    /// return interpolation model
    virtual QualityType fit1d(const RawDataArrayType& /* range */, std::unique_ptr<InterpolationModel>& /* model */);



  protected:
    /// standard derivation in bounding box
    CoordinateType tolerance_stdev_box_;
    /// basic statistics
    Math::BasicStatistics<> statistics_;
    /// interpolation step size
    CoordinateType interpolation_step_;

    void updateMembers_() override;
  };

} // namespace OpenMS
