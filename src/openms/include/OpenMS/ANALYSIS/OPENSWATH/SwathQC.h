// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/OpenMSConfig.h>

#include <functional>
#include <map>
#include <vector>

namespace OpenMS
{
  class String;
  class MSSpectrum;
  class ExperimentalSettings;

}

namespace OpenSwath
{

  /**
    @brief Quality Control function for OpenSwath

    The class is meant to collect information as a plugin of type IMSDataConsumer.
    You can (if the relevant data is already available and you do not mind the overhead), call
    static functions, which will do the same.

  */
  class OPENMS_DLLAPI SwathQC
  {
  public:
    typedef std::map<int,int> ChargeDistribution;
    
    /// default C'tor (forbidden)
    SwathQC() = delete;

    /**
      @brief CTor with arguments
      
      @param cd_spectra Number of MS spectra to inspect for charge distribution estimation
      @param decon_ms1_mz_tol m/z tolerance for isotope deconvolution
    */
    SwathQC(const size_t cd_spectra, const double decon_ms1_mz_tol);


    /**
      @brief Returns a lambda function, which captures internal members and can be used in MSDataTransformingConsumer::setExperimentalSettingsFunc()


       If @p es contains a metavalue 'nr_ms1_spectra' it will be used to set the internal number of expected MS1 spectra
    */
    std::function<void(const OpenMS::ExperimentalSettings& es)> getExpSettingsFunc();

    /// returns a lambda function, which captures internal members and can be used in MSDataTransformingConsumer::setExperimentalSettingsFunc()
    /// PeakPicking is performed internally if the data is estimated to be profile data.
    /// Uses the Deisotoper class for charge determination and updates the internal charge count.
    std::function<void (const OpenMS::MSSpectrum&)> getSpectraProcessingFunc();


    /**
      @brief Sample the spectra in all MS1 Swath maps and determine all charge states.

      Conveniently wraps internal functions and applies them to a whole set of data.
      However, it might be more speed efficient to sample the spectra as they are loaded using member functions.

      From all swath maps in @p swath_maps which are ms_level == 1 we extract @p nr_samples spectra (subsampling),
      determine the charge states of all isotopic envelopes and return their total counts
      using getSpectraProcessingFunc().

      @param[in] swath_maps Swath maps of mixed ms-level
      @param[in] nr_samples Number of spectra to sample per Swath map. To sample all spectra, set to -1
      @param[in] mz_tol     Error tolerance in Th in which an isotopic peak is expected (assuming C12-C13 distance)
      @return Distribution of charge (key = charge, value = counts)

      @throw Exception::Postcondition if Deisotoper did not return charge data

    */
    static ChargeDistribution getChargeDistribution(const std::vector<SwathMap>& swath_maps, const size_t nr_samples, const double mz_tol);


    /**
      @brief Save all internally collected data to a JSON file.
    */
    void storeJSON(const OpenMS::String& filename);

    /**
      @brief returns the charge distribution which was internally computed by applying getSpectraProcessingFunc() externally
    */
    const ChargeDistribution& getChargeDistribution() const;

    /**
      @brief Explicitly set the number of expected MS1 spectra (for sampling charge distribution)

      Computing the charge distribution using getSpectraProcessingFunc() requires knowing the total number of MS1 spectra.
      Either use getExpSettingsFunc() externally, or use this method to set it explicitly (depending on workflow).
      If @p nr is set to 0, all spectra passed into getSpectraProcessingFunc() will be inspected for their charge distribution.
    */
    void setNrMS1Spectra(size_t nr);


   protected: 

    /**
      @brief Given a total spectrum count and a number of spectra to inspect, is the current index a candidate?

      Allows to uniformly sample a range of spectra.
      E.g. given 10 spectra exist, and we want to subsample 4 of them, the answer of this function for every
      requested index from 0..9 is:

      If the total number of spectra is unknown, pass 0 as first argument, which will return true for every query, i.e. sample everything.

      @param[in] total_spec_count Total number of spectra expected
      @param[in] subsample_count Number of spectra which should be sampled from @p total_spec_count
      @param[in] idx Index of the spectrum under question

    */
    static bool isSubsampledSpectrum_(const size_t total_spec_count, const size_t subsample_count, const size_t idx);
   
   private: 

    /// internal ChargeDistribution which is augmented upon calling the corresponding member functions
    ChargeDistribution cd_;

    /// number of MS1 spectra expected
    size_t nr_ms1_spectra_;
    /// number of spectra to inspect for charge distribution
    size_t cd_spectra_;
    /// m/z tolerance for isotope deconvolution
    double decon_ms1_mz_tol_;

    /// keeps track of number of spectra passed to getSpectraProcessingFunc()
    size_t ms1_spectra_seen_;
  };
    
}

