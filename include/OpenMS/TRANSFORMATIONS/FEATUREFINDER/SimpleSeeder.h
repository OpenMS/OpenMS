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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include <algorithm>
#include <vector>
#include <iostream>

namespace OpenMS
{
  /**
      @brief Simple seeding class that uses the strongest peak as next seed.

      This class simply sorts the peaks according to intensity and proposes
      the highest peak, which is not yet included in a feature, as next seed.

      @htmlinclude OpenMS_SimpleSeeder.parameters

      @ingroup FeatureFinder
  */
  template <class PeakType, class FeatureType>
  class SimpleSeeder :
    public FeaFiModule<PeakType, FeatureType>,
    public FeatureFinderDefs
  {
public:
    typedef FeaFiModule<PeakType, FeatureType> Base;
    typedef MSExperiment<PeakType> MapType;

    /// Constructor
    SimpleSeeder(const MSExperiment<PeakType> * map, FeatureMap<FeatureType> * features, FeatureFinder * ff) :
      Base(map, features, ff),
      initialized_(false)
    {
      this->setName("SimpleSeeder");

      this->defaults_.setValue("min_intensity", 0.0, "Absolute value for the minimum intensity required for a seed.");
      this->defaults_.setMinFloat("min_intensity", 0.0);
      this->defaults_.setValue("signal_to_noise", 10.0, "Minimal required SignalToNoise (S/N) ratio for a seed.");
      this->defaults_.setMinFloat("signal_to_noise", 0.0);

      //this->subsections_.push_back("SignalToNoiseEstimationParameter");
      SignalToNoiseEstimatorMedian<typename MapType::SpectrumType> sne;              // make sure this is the same as in pick()!
      this->defaults_.insert("SignalToNoiseEstimationParameter:", sne.getDefaults());

      this->defaultsToParam_();
    }

    /// destructor
    virtual ~SimpleSeeder()
    {
    }

    /// return the next seed
    IndexPair nextSeed()
    {
      if (!initialized_)
      {
        initialize_();
      }

      // while the current peak is either already used or in a feature jump to next peak...
      while (current_peak_ != indices_.end() && this->ff_->getPeakFlag(*current_peak_) == USED)
      {
        ++current_peak_;
      }

      if (current_peak_ == indices_.end())
      {
        // if no seed was found:
        if (indices_.empty()) throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, IndexPair());
        else throw NoSuccessor(__FILE__, __LINE__, __PRETTY_FUNCTION__, *(current_peak_ - 1));
      }

      this->ff_->setProgress(current_peak_ - indices_.begin());

      // set flag
      this->ff_->getPeakFlag(*current_peak_) = USED;

      return *(current_peak_++);
    }         // nextSeed

protected:

    void initialize_()
    {
      // determine minimum intensity and signal-to-noise parameter for last seed
      typename FeatureType::IntensityType noise_threshold  = this->param_.getValue("min_intensity");
      typename FeatureType::IntensityType sn  = this->param_.getValue("signal_to_noise");

#ifdef DEBUG_FEATUREFINDER
      std::cout << "Intensity threshold: " << noise_threshold << std::endl;
      std::cout << "S/N: " << sn << std::endl;
#endif

      // fill indices_ for peaks above noise threshold and S/N
      IndexPair tmp = std::make_pair(0, 0);
      if (sn == 0)
      {
        while (tmp.first < (*this->map_).size())
        {
          tmp.second = 0;
          while (tmp.second < (*this->map_)[tmp.first].size())
          {
            if (this->getPeakIntensity(tmp) > noise_threshold)
            {
              indices_.push_back(tmp);
            }
            ++tmp.second;
          }
          ++tmp.first;
        }
      }
      else
      {
        SignalToNoiseEstimatorMedian<typename MapType::SpectrumType> estimator;
        Param param(this->param_.copy("SignalToNoiseEstimationParameter:", true));
        estimator.setParameters(param);

        for (typename MapType::ConstIterator it = (*this->map_).begin(); it != (*this->map_).end(); ++it)
        {
          estimator.init(it->begin(), it->end());
          tmp.second = 0;
          for (typename MapType::SpectrumType::ConstIterator spec = it->begin(); spec != it->end(); ++spec)
          {
            if (estimator.getSignalToNoise(spec) > sn && this->getPeakIntensity(tmp) > noise_threshold)
            {
              indices_.push_back(tmp);
            }
            ++tmp.second;
          }
          ++tmp.first;
        }
      }

#ifdef DEBUG_FEATUREFINDER
      std::cout << "Number of peaks above threshold (" << noise_threshold   << ") and S/N (" << sn << "): " << indices_.size() << std::endl;
#endif

      // sort index vector by intensity of peaks (highest first)
      sort(indices_.begin(), indices_.end(),
           reverseComparator(Internal::IntensityLess<Base>(*this))
           );

      // progress logger
      this->ff_->startProgress(0, indices_.size(), "FeatureFinder");

      current_peak_ = indices_.begin();

      initialized_ = true;
    }

    /// contains the indices
    std::vector<IndexPair> indices_;

    /// Points to the next peak in the peak vector
    std::vector<IndexPair>::const_iterator current_peak_;

    /// Flag that indicates of the indices are initialized
    bool initialized_;

private:
    /// Not implemented
    SimpleSeeder();
    /// Not implemented
    SimpleSeeder & operator=(const SimpleSeeder &);
    /// Not implemented
    SimpleSeeder(const SimpleSeeder &);

  };   // class SimpleSeeder

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
