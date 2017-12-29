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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  // forward declaration
  class FeatureFinder;
  class FeatureMap;

  /// Summary of fitting results
  struct OPENMS_DLLAPI Summary
  {
    std::map<String, UInt> exception; //count exceptions
    UInt no_exceptions;
    std::map<String, UInt> mz_model; //count used mz models
    std::map<float, UInt> mz_stdev; //count used mz standard deviations
    std::vector<UInt> charge; //count used charges
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
    public DefaultParamHandler
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
      features_(nullptr),
      ff_(nullptr)
    {
    }

    /// destructor
    ~FeatureFinderAlgorithm() override
    {
    }

    /// register all derived classes here (see FeatureFinderAlgorithm_impl.h)
    static void registerChildren();

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
    void setData(const MapType& map, FeatureMap& features, FeatureFinder& ff)
    {
      map_ = &map;
      features_ = &features;
      ff_ = &ff;
    }

    /**
        @brief Sets a reference to the calling FeatureFinder

        @exception Exception::IllegalArgument is thrown if the algorithm does not support user-specified seed lists
    */
    virtual void setSeeds(const FeatureMap& seeds)
    {
      if (seeds.size() != 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The used feature detection algorithm does not support user-specified seed lists!");
      }
    }

protected:

    /// Input data pointer
    const MapType* map_;

    /// Output data pointer
    FeatureMap* features_;

    /// Pointer to the calling FeatureFinder that is used to access the feature flags
    FeatureFinder* ff_;

private:

    /// Not implemented
    FeatureFinderAlgorithm& operator=(const FeatureFinderAlgorithm&);

    /// Not implemented
    FeatureFinderAlgorithm(const FeatureFinderAlgorithm&);

  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
