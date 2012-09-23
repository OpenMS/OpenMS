// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>

namespace OpenMS
{
  /**
  @brief FeatureFinderAlgorithm implementation using the Simple* modules.

  @deprecated Deprecated in OpenMS 1.7.

  SimpleSeeder, SimpleExtender, ModelFitter (using EmgModel (exponentially
  modified gaussian with parameter optimization [Levenberg-Marquardt
  algorithm]) in RT dimension and IsotopeModel (charge does not equal zero)
  or LmaGaussModel (parameter optimization using Levenberg-Marquardt
  algorithm) in dimension of mz).

  @htmlinclude OpenMS_FeatureFinderAlgorithmSimple.parameters

  @ingroup FeatureFinder
  */
  template <class PeakType, class FeatureType>
  class FeatureFinderAlgorithmSimple :
    public FeatureFinderAlgorithm<PeakType, FeatureType>,
    public FeatureFinderDefs
  {

public:
    /// default constructor
    FeatureFinderAlgorithmSimple() :
      FeatureFinderAlgorithm<PeakType, FeatureType>()
    {
      this->defaults_ = getDefaultParameters();
      this->check_defaults_ =  false;
    }

    virtual Param getDefaultParameters() const
    {
      Param tmp;

      SimpleSeeder<PeakType, FeatureType> seeder(this->map_, this->features_, this->ff_);
      tmp.insert("seeder:", seeder.getParameters());
      tmp.setSectionDescription("seeder", "Settings for the seeder (Determines potential feature regions)");

      SimpleExtender<PeakType, FeatureType> extender(this->map_, this->features_, this->ff_);
      tmp.insert("extender:", extender.getParameters());
      tmp.setSectionDescription("extender", "Settings for the extender (Collects all peaks belonging to a feature)");

      ModelFitter<PeakType, FeatureType> fitter(this->map_, this->features_, this->ff_);
      tmp.insert("fitter:", fitter.getParameters());
      tmp.setSectionDescription("fitter", "Settings for the modefitter (Fits a model to the data determinging the probapility that they represent a feature.)");

      return tmp;
    }

    virtual void run()
    {
#ifdef DEBUG_FEATUREFINDER
      UInt seed_nr = 0;
#endif
      SimpleSeeder<PeakType, FeatureType> seeder(this->map_, this->features_, this->ff_);
      seeder.setParameters(this->getParameters().copy("seeder:", true));

      SimpleExtender<PeakType, FeatureType> extender(this->map_, this->features_, this->ff_);
      extender.setParameters(this->getParameters().copy("extender:", true));

      ModelFitter<PeakType, FeatureType> fitter(this->map_, this->features_, this->ff_);
      Param params;
      params.setDefaults(this->getParameters().copy("fitter:", true));
      params.setValue("fit_algorithm", "simple");
      fitter.setParameters(params);

      /// Summary of fitting results
      Summary summary;

      try
      {
        for (;; )
        {
#ifdef DEBUG_FEATUREFINDER
          std::cout << "===============================" << std::endl;
          std::cout << "### Seeder (seed # " << ++seed_nr << ")..." << std::endl;
#endif
          IndexPair seed = seeder.nextSeed();

#ifdef DEBUG_FEATUREFINDER
          std::cout << "seed ... " << seed.first << " - " << seed.second << std::endl;
          std::cout << "### Extender..." << std::endl;
#endif
          ChargedIndexSet index_set;
          index_set.insert(seed);
          ChargedIndexSet region;
          extender.extend(index_set, region);

#ifdef DEBUG_FEATUREFINDER
          std::cout << "### ModelFitter..." << std::endl;
#endif
          try
          {
            this->features_->push_back(fitter.fit(region));

            // gather information for fitting summary
            {
              const Feature & f = this->features_->back();

              // quality, correlation
              DoubleReal corr = f.getOverallQuality();
              summary.corr_mean += corr;
              if (corr < summary.corr_min) summary.corr_min = corr;
              if (corr > summary.corr_max) summary.corr_max = corr;

              // charge
              UInt ch = f.getCharge();
              if (ch >= summary.charge.size())
              {
                summary.charge.resize(ch + 1);
              }
              summary.charge[ch]++;

              // MZ model type
              const Param & p = f.getModelDescription().getParam();
              ++summary.mz_model[p.getValue("MZ")];

              // standard deviation of isotopic peaks
              if (p.exists("MZ:isotope:stdev") && p.getValue("MZ:isotope:stdev") != DataValue::EMPTY)
              {
                ++summary.mz_stdev[p.getValue("MZ:isotope:stdev")];
              }
            }
          }
          catch (Exception::UnableToFit ex)
          {
            //std::cout << "UnableToFit: " << ex.what() << std::endl;

            // set unused flag for all data points
            for (IndexSet::const_iterator it = region.begin(); it != region.end(); ++it)
            {
              this->ff_->getPeakFlag(*it) = UNUSED;
            }

            // gather information for fitting summary
            {
              ++summary.no_exceptions;
              ++summary.exception[ex.getName()];
            }
          }
        }         // for
      }       // try
      catch (NoSuccessor ex)
      {
      }

      this->ff_->endProgress();

      // print fitting summary
      {
        Size size = this->features_->size();
        std::cout << size << " features were found. " << std::endl;

        // compute corr_mean
        summary.corr_mean /= size;

        std::cout << "FeatureFinder summary:\n"
                  << "Correlation:\n\tminimum: " << summary.corr_min << "\n\tmean: " << summary.corr_mean
                  << "\n\tmaximum: " << summary.corr_max << std::endl;

        std::cout << "Exceptions:\n";
        for (std::map<String, UInt>::const_iterator it = summary.exception.begin(); it != summary.exception.end(); ++it)
        {
          std::cout << "\t" << it->first << ": " << it->second * 100 / summary.no_exceptions << "% (" << it->second << ")\n";
        }

        std::cout << "Chosen mz models:\n";
        for (std::map<String, UInt>::const_iterator it = summary.mz_model.begin(); it != summary.mz_model.end(); ++it)
        {
          std::cout << "\t" << it->first << ": " << it->second * 100 / size << "% (" << it->second << ")\n";
        }

        std::cout << "Chosen mz stdevs:\n";
        for (std::map<float, UInt>::const_iterator it = summary.mz_stdev.begin(); it != summary.mz_stdev.end(); ++it)
        {
          std::cout << "\t" << it->first << ": " << it->second * 100 / (size - summary.charge[0]) << "% (" << it->second << ")\n";
        }

        std::cout << "Charges:\n";
        for (Size i = 1; i < summary.charge.size(); ++i)
        {
          if (summary.charge[i] != 0)
          {
            std::cout << "\t+" << i << ": " << summary.charge[i] * 100 / (size - summary.charge[0]) << "% (" << summary.charge[i] << ")\n";
          }
        }
      }
    }     // run

    static FeatureFinderAlgorithm<PeakType, FeatureType> * create()
    {
      return new FeatureFinderAlgorithmSimple();
    }

    static const String getProductName()
    {
      return "simple";
    }

private:
    /// Not implemented
    FeatureFinderAlgorithmSimple & operator=(const FeatureFinderAlgorithmSimple &);
    /// Not implemented
    FeatureFinderAlgorithmSimple(const FeatureFinderAlgorithmSimple &);

  };   // FeatureFinderAlgorithmSimple

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSIMPLE_H
