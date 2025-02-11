// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Witold Wolski, Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <utility>

#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS::DIAHelpers
  {

    void adjustExtractionWindow(double& right, double& left, const double& mz_extract_window, const bool& mz_extraction_ppm)
    {
      OPENMS_PRECONDITION(mz_extract_window > 0, "MZ extraction window needs to be larger than zero.");

      if (mz_extraction_ppm)
      {
        left -= left * mz_extract_window / 2e6;
        right += right * mz_extract_window / 2e6;
      }
      else
      {
        left -= mz_extract_window / 2.0;
        right += mz_extract_window / 2.0;
      }
      // If the left value is now < 0, this is invalid correct it to be 0
      if (left < 0)
      {
        left = 0;
      }
    }

    // Helper for integrate window returns the sum of all intensities, sum of all ion mobilities and sum of all mz
    // no expensive division calls
    // assumes mz, im and intensity should already be initiated.
    void integrateWindow_(const OpenSwath::SpectrumPtr& spectrum,
                double & mz,
                double & im,
                double & intensity,
                const RangeMZ & mz_range,
                const RangeMobility & im_range,
                bool centroided)
    {
      OPENMS_PRECONDITION(spectrum != nullptr, "precondition: Spectrum cannot be nullptr");
      OPENMS_PRECONDITION(spectrum->getMZArray() != nullptr, "precondition: Cannot integrate if no m/z is available.");
      //OPENMS_PRECONDITION(!spectrum->getMZArray()->data.empty(), " precondition: Warning: Cannot integrate if spectrum is empty"); // This is not a failure should check for this afterwards
      OPENMS_PRECONDITION(std::adjacent_find(spectrum->getMZArray()->data.begin(),
        spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
        "Precondition violated: m/z vector needs to be sorted!" );
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "precondition: MZ and Intensity array need to have the same length.");

      // ion mobility specific preconditions
      //OPENMS_PRECONDITION((im_range.isEmpty()) && (spectrum->getDriftTimeArray() != nullptr), "precondition: Cannot integrate with drift time if no drift time is available."); This is not a failure can handle this
      OPENMS_PRECONDITION((spectrum->getDriftTimeArray() == nullptr) || (spectrum->getDriftTimeArray()->data.empty()) || (spectrum->getMZArray()->data.size() == spectrum->getDriftTimeArray()->data.size()), "precondition: MZ and Drift Time array need to have the same length.");
      OPENMS_PRECONDITION(!centroided,  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION));

      if ( spectrum->getMZArray()->data.empty() )
      {
        OPENMS_LOG_WARN << "Warning: Cannot integrate if spectrum is empty" << std::endl;
        return;
      }

      // if im_range is set, than integrate across dirft time
      if (!im_range.isEmpty()) // if imRange supplied, integrate across IM
      {
        if (spectrum->getDriftTimeArray() == nullptr)
        {
            //throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot integrate with drift time if no drift time is available");
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot integrate with drift time if no drift time is available");
        }
      }

      if (!centroided)
      {
        // get the weighted average for noncentroided data.
        // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
        auto mz_arr_end = spectrum->getMZArray()->data.end();
        auto int_it = spectrum->getIntensityArray()->data.begin();

        // this assumes that the spectra are sorted!
        auto mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(), mz_arr_end, mz_range.getMin());

        // also advance intensity iterator now
        auto iterator_pos = std::distance(spectrum->getMZArray()->data.begin(), mz_it);
        std::advance(int_it, iterator_pos);

        double mz_end = mz_range.getMax(); // store the maximum mz value in a double to minimize function calls

        if ( !im_range.isEmpty() ) // integrate across im as well
        {
          auto im_it = spectrum->getDriftTimeArray()->data.begin();

          // also advance ion mobility iterator now
          std::advance(im_it, iterator_pos);

          // Start iteration from mz start, end iteration when mz value is larger than mz_end, only store only storing ion mobility values that are in the range
          //while ( (mz_it != mz_arr_end) && (*mz_it < mz_end) )
          while ( (mz_it != mz_arr_end) && (*mz_it < mz_end) )
          {
            if (im_range.contains(*im_it))
            {
              intensity += (*int_it);
              im += (*int_it) * (*im_it);
              mz += (*int_it) * (*mz_it);
            }
            ++mz_it;
            ++int_it;
            ++im_it;
          }
        }
        else // where do not have IM
        {
          while ( mz_it != mz_arr_end && *mz_it < mz_end )
          {
            intensity += (*int_it);
            mz += (*int_it) * (*mz_it);

            ++mz_it;
            ++int_it;
          }
        }
      }
      else
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }


    void integrateWindows(const OpenSwath::SpectrumPtr& spectrum,
                          const std::vector<double> & windowsCenter,
                          double width,
                          std::vector<double> & integratedWindowsIntensity,
                          std::vector<double> & integratedWindowsMZ,
			  std::vector<double> & integratedWindowsIm,
                          const RangeMobility & range_im,
                          bool remZero)
    {
      std::vector<double>::const_iterator beg = windowsCenter.begin();
      std::vector<double>::const_iterator end = windowsCenter.end();
      double mz, intensity, im;
      for (; beg != end; ++beg)
      {
        // assemble RangeMZ object based on window
        RangeMZ range_mz(*beg);
        range_mz.minSpanIfSingular(width);

        if (integrateWindow(spectrum, mz, im, intensity, range_mz, range_im, false))
        {
          integratedWindowsIntensity.push_back(intensity);
          integratedWindowsMZ.push_back(mz);
	  integratedWindowsIm.push_back(im);
        }
        else if (!remZero)
        {
          integratedWindowsIntensity.push_back(0.);
          integratedWindowsMZ.push_back(*beg);
          if ( !range_im.isEmpty() )
          {
            integratedWindowsIm.push_back( range_im.center() ); // average drift time
          }
          else
          {
            integratedWindowsIm.push_back(-1);
          }
        }
      }
    }


    void integrateWindows(const SpectrumSequence& spectra,
                          const std::vector<double> & windowsCenter,
                          double width,
                          std::vector<double> & integratedWindowsIntensity,
                          std::vector<double> & integratedWindowsMZ,
			  std::vector<double> & integratedWindowsIm,
                          const RangeMobility& range_im,
                          bool remZero)
    {

      double mz(-1), intensity(0), im(-1);
      if (windowsCenter.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No windows supplied!");
        return;
      }

      if (spectra.empty())
      {
        OPENMS_LOG_WARN << "Warning: no spectra provided" << std::endl;
        return;
      }

      std::vector<double>::const_iterator beg = windowsCenter.begin();
      std::vector<double>::const_iterator end = windowsCenter.end();
      for (; beg != end; ++beg)
      {
        //assemble rangeMZ object based on windows
        RangeMZ range_mz(*beg);
        range_mz.minSpanIfSingular(width);

        if (integrateWindow(spectra, mz, im, intensity, range_mz, range_im, false))
        {
          integratedWindowsIntensity.push_back(intensity);
          integratedWindowsMZ.push_back(mz);
	  integratedWindowsIm.push_back(im);
        }
        else if (!remZero)
        {
          integratedWindowsIntensity.push_back(0.);
          integratedWindowsMZ.push_back(*beg); // push back center of window
          if ( !range_im.isEmpty() )
          {
            integratedWindowsIm.push_back(range_im.center()); // push back average drift
          }
          else
          {
            integratedWindowsIm.push_back(-1);
          }
        }
      }
    }

    bool integrateWindow(const OpenSwath::SpectrumPtr& spectrum,
                                              double & mz,
                                              double & im,
                                              double & intensity,
                                              const RangeMZ & range_mz,
                                              const RangeMobility & range_im,
                                              bool centroided)
    {

      // initiate the values
      mz = 0;
      im = 0;
      intensity = 0;

      integrateWindow_(spectrum, mz, im, intensity, range_mz, range_im, centroided);

      // Post processing get the weighted average mz and im by dividing my intensity
      if (intensity > 0.)
      {
        mz /= intensity;

        if ( !range_im.isEmpty() )
        {
          im /= intensity;
        }
        else
        {
          im = -1;
        }
        return true;
      }
      else
      {
        im = -1;
        mz = -1;
        intensity = 0;
        return false;
      }
    }

    bool integrateWindow(const SpectrumSequence& spectra,
                                              double & mz,
                                              double & im,
                                              double & intensity,
                                              const RangeMZ & range_mz,
                                              const RangeMobility & range_im,
                                              bool centroided)
    {
      // initiate the values
      mz = 0;
      im = 0;
      intensity = 0;

      if (!spectra.empty())
      {
        for (const auto& s : spectra)
        {
          integrateWindow_(s, mz, im, intensity, range_mz, range_im, centroided);
        }

        // Post processing get the weighted average mz and im by dividing my intensity
        if (intensity > 0.)
        {
          mz /= intensity;
          if ( !range_im.isEmpty() )
          {
            im /= intensity;
          }
          else // if no IM set to -1
          {
            im = -1;
          }
          return true;
        }
        else
        {
          // if (intensity <= 0)
          im = -1;
          mz = -1;
          intensity = 0;
          return false;
        }
      }
      else
      {
        // if (all_spectra.empty())
        OPENMS_LOG_WARN << "Warning: no spectra provided" << std::endl;
        im = -1;
        mz = -1;
        intensity = 0;
        return false;
      }
    }

    // for SWATH -- get the theoretical b and y series masses for a sequence
    void getBYSeries(const AASequence& a, //
                     std::vector<double>& bseries, //
                     std::vector<double>& yseries, //
                     TheoreticalSpectrumGenerator const * generator,
                     int charge)
    {
      // Note: We pass TheoreticalSpectrumGenerator ptr, as constructing it each time is too slow.
      OPENMS_PRECONDITION(charge > 0, "For constructing b/y series we require charge being a positive integer");

      if (a.empty()) return;

      PeakSpectrum spec;
      generator->getSpectrum(spec, a, charge, charge);

      // Data array is present if AASequence is not empty
      const PeakSpectrum::StringDataArray& ion_name = spec.getStringDataArrays()[0];

      for (Size i = 0; i != spec.size(); ++i)
      {
        if (ion_name[i][0] == 'y')
        {
          yseries.push_back(spec[i].getMZ());
        }
        else if (ion_name[i][0] == 'b')
        {
          bseries.push_back(spec[i].getMZ());
        }
      }
    } // end getBYSeries

    // for SWATH -- get the theoretical b and y series masses for a sequence
    void getTheorMasses(const AASequence& a,
                        std::vector<double>& masses,
                        TheoreticalSpectrumGenerator const * generator,
                        int charge)
    {
      // Note: We pass TheoreticalSpectrumGenerator ptr, as constructing it each time is too slow.
      OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer");

      PeakSpectrum spec;
      generator->getSpectrum(spec, a, charge, charge);
      for (PeakSpectrum::iterator it = spec.begin();
           it != spec.end(); ++it)
      {
        masses.push_back(it->getMZ());
      }
    } // end getBYSeries

    void  getAveragineIsotopeDistribution(const double product_mz,
                                         std::vector<std::pair<double, double> >& isotopes_spec,
                                         int charge,
                                         const int nr_isotopes,
                                         const double mannmass)
    {
      charge = std::abs(charge);
      typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;
      // create the theoretical distribution
      CoarseIsotopePatternGenerator solver(nr_isotopes);
      TheoreticalIsotopePattern isotopes;
      //Note: this is a rough estimate of the weight, usually the protons should be deducted first, left for backwards compatibility.
      auto d = solver.estimateFromPeptideWeight(product_mz * charge);

      double mass = product_mz;
      for (IsotopeDistribution::Iterator it = d.begin(); it != d.end(); ++it)
      {
        isotopes_spec.emplace_back(mass, it->getIntensity());
        mass += mannmass / charge;
      }
    } //end of dia_isotope_corr_sub

    //simulate spectrum from AASequence
    void simulateSpectrumFromAASequence(const AASequence& aa,
                                        std::vector<double>& first_isotope_masses, //[out]
                                        std::vector<std::pair<double, double> >& isotope_masses, //[out]
                                        TheoreticalSpectrumGenerator const * generator, int charge)
    {
      getTheorMasses(aa, first_isotope_masses, generator, charge);
      for (std::size_t i = 0; i < first_isotope_masses.size(); ++i)
      {
        getAveragineIsotopeDistribution(first_isotope_masses[i], isotope_masses,
                                        charge);
      }
    }

    /// given an experimental spectrum add isotope pattern.
    void addIsotopes2Spec(const std::vector<std::pair<double, double> >& spec,
                          std::vector<std::pair<double, double> >& isotope_masses, //[out]
                          Size nr_isotopes, int charge)
    {

      for (std::size_t i = 0; i < spec.size(); ++i)
      {
        std::vector<std::pair<double, double> > isotopes;
        getAveragineIsotopeDistribution(spec[i].first, isotopes, charge, nr_isotopes);
        for (Size j = 0; j < isotopes.size(); ++j)
        {
          isotopes[j].second *= spec[i].second; //multiple isotope intensity by spec intensity
          isotope_masses.push_back(isotopes[j]);
        }
      }
    }

    /// given a peak of experimental mz and intensity, add isotope pattern to a "spectrum".
    void addSinglePeakIsotopes2Spec(double mz, double ity,
                                    std::vector<std::pair<double, double> >& isotope_masses, //[out]
                                    Size nr_isotopes, int charge)
    {
      std::vector<std::pair<double, double> > isotopes;
      getAveragineIsotopeDistribution(mz, isotopes, charge, nr_isotopes);
      for (Size j = 0; j < isotopes.size(); ++j)
      {
        isotopes[j].second *= ity; //multiple isotope intensity by spec intensity
        isotope_masses.push_back(isotopes[j]);
      }
    }

    //Add masses before first isotope
    void addPreisotopeWeights(const std::vector<double>& first_isotope_masses,
                              std::vector<std::pair<double, double> >& isotope_spec, // output
                              UInt nr_peaks, double pre_isotope_peaks_weight, // weight of pre isotope peaks
                              double mannmass, int charge)
    {
      charge = std::abs(charge);
      for (std::size_t i = 0; i < first_isotope_masses.size(); ++i)
      {
        Size mul = 1.;
        for (UInt j = 0; j < nr_peaks; ++j, ++mul)
        {
          isotope_spec.emplace_back(first_isotope_masses[i] - (mul * mannmass) / charge,
                                    pre_isotope_peaks_weight);
        }
      }
      sortByFirst(isotope_spec);
    }

    //Add masses before first isotope
    void addPreisotopeWeights(double mz,
                              std::vector<std::pair<double, double> >& isotope_spec, // output
                              UInt nr_peaks, double pre_isotope_peaks_weight, // weight of pre isotope peaks
                              double mannmass, int charge)
    {
      charge = std::abs(charge);
      Size mul = 1.;
      for (UInt j = 0; j < nr_peaks; ++j, ++mul)
      {
        isotope_spec.emplace_back(mz - (mul * mannmass) / charge,
                                  pre_isotope_peaks_weight);
      }
    }

    struct MassSorter
    {
      bool operator()(const std::pair<double, double>& left,
                      const std::pair<double, double>& right)
      {
        return left.first < right.first;
      }
    };

    void sortByFirst(std::vector<std::pair<double, double> >& tmp)
    {
      std::sort(tmp.begin(), tmp.end(), MassSorter());
    }

    void extractFirst(const std::vector<std::pair<double, double> >& peaks,
                      std::vector<double>& mass)
    {
      for (std::size_t i = 0; i < peaks.size(); ++i)
      {
        mass.push_back(peaks[i].first);
      }
    }

    void extractSecond(const std::vector<std::pair<double, double> >& peaks,
                       std::vector<double>& mass)
    {
      for (std::size_t i = 0; i < peaks.size(); ++i)
      {
        mass.push_back(peaks[i].second);
      }
    }

    RangeMZ createMZRangePPM(const double mzRef, const double dia_extraction_window, const bool is_ppm)
    {
      RangeMZ rangeMZ(mzRef);
      if (is_ppm)
      {
        double ppm =  Math::ppmToMass(dia_extraction_window, mzRef);
        rangeMZ.minSpanIfSingular(ppm);
      }
      else
      {
        rangeMZ.minSpanIfSingular(dia_extraction_window);
      }
      return rangeMZ;
    }

  }
