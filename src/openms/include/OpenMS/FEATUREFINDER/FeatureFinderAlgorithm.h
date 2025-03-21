// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


namespace OpenMS
{

  // forward declaration
  class FeatureMap;

  /// Summary of fitting results
  struct OPENMS_DLLAPI Summary
  {
    std::map<String, UInt> exception; ///<count exceptions
    UInt no_exceptions;
    std::map<String, UInt> mz_model; ///<count used mz models
    std::map<float, UInt> mz_stdev; ///<count used mz standard deviations
    std::vector<UInt> charge; ///<count used charges
    double corr_mean, corr_max, corr_min; //boxplot for correlation

    /// Initial values
    Summary() :
      no_exceptions(0),
      corr_mean(0),
      corr_max(0),
      corr_min(1)
    {}

  };

  /**
      @brief Abstract base class for FeatureFinder algorithms

  */
  class FeatureFinderAlgorithm :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Input map type
    typedef PeakMap MapType;
    /// Coordinate/Position type of peaks
    typedef MapType::CoordinateType CoordinateType;
    /// Intensity type of peaks
    typedef MapType::IntensityType IntensityType;

    /// default constructor
    FeatureFinderAlgorithm() :
      DefaultParamHandler("FeatureFinderAlgorithm"),
      map_(nullptr),
      features_(nullptr)
    {
    }

    /// destructor
    ~FeatureFinderAlgorithm() override
    {
    }

    /// register all derived classes here (see FeatureFinderAlgorithm_impl.h)


    /// Main method that implements the actual algorithm
    virtual void run() = 0;

    /**
        @brief Returns the default parameters. Reimplement

        Reimplement if you derive a class and have to incorporate sub-algorithm default parameters.
    */
    virtual Param getDefaultParameters() const
    {
      return this->defaults_;
    }

    /// Sets a reference to the calling FeatureFinder
    void setData(const MapType& map, FeatureMap& features)
    {
      map_ = &map;
      features_ = &features;
    }

    /**
        @brief Sets a reference to the calling FeatureFinder

        @exception Exception::IllegalArgument is thrown if the algorithm does not support user-specified seed lists
    */
    virtual void setSeeds(const FeatureMap& seeds)
    {
      if (!seeds.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The used feature detection algorithm does not support user-specified seed lists!");
      }
    }

protected:

    /// Input data pointer
    const MapType* map_;

    /// Output data pointer
    FeatureMap* features_;

private:

    /// Not implemented
    FeatureFinderAlgorithm& operator=(const FeatureFinderAlgorithm&);

    /// Not implemented
    FeatureFinderAlgorithm(const FeatureFinderAlgorithm&);

  };
}

