// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  class ConsensusMap;
  class ConsensusFeature;

  /**
    @brief Extracts individual channels from MS/MS spectra for isobaric labeling experiments.

    In addition to extracting the channel information this class can also filter the extracted channel
    information according to several parameters, i.e., discard channel information if certain criteria
    are not met.

    @li %Precursor activation method (e.g., select only HCD scans).
    @li Minimum precursor intensity.
    @li Minimum reporter intensity (i.e., remove reporter channels below a certain intensity)
    @li %Precursor purity (i.e., fraction of TIC in the precursor window that can be assigned to the precursor)

    The precursor purity computation uses the interpolation approach described in:
    Savitski MM, Sweetman G, Askenazi M, et al. (2011). Delayed fragmentation and optimized isolation width settings
    for improvement of protein identification and accuracy of isobaric mass tag quantification on Orbitrap-type mass
    spectrometers. Analytical chemistry 83: 8959-67. http://www.ncbi.nlm.nih.gov/pubmed/22017476

    @note Centroided MS and MS/MS data is required.

    @htmlinclude OpenMS_IsobaricChannelExtractor.parameters
  */
  class OPENMS_DLLAPI IsobaricChannelExtractor :
    public DefaultParamHandler
  {
public:
    /**
      @brief C'tor to create a new channel extractor for the given quantitation method.

      @param quant_method IsobaricQuantitationMethod providing the necessary information which channels should be extracted.
    */
    explicit IsobaricChannelExtractor(const IsobaricQuantitationMethod* const quant_method);

    /// Copy c'tor
    IsobaricChannelExtractor(const IsobaricChannelExtractor& other);

    /// Assignment operator
    IsobaricChannelExtractor& operator=(const IsobaricChannelExtractor& rhs);

    /**
      @brief Extracts the isobaric channels from the tandem MS data and stores intensity values in a consensus map.

      @param ms_exp_data Raw data to search for isobaric quantitation channels.
      @param consensus_map Output map containing the identified channels and the corresponding intensities.
    */
    void extractChannels(const PeakMap& ms_exp_data, ConsensusMap& consensus_map);

private:
    /**
      @brief Small struct to capture the current state of the purity computation.

      It basically contains two iterators pointing to the current potential
      MS1 precursor scan of an MS2 scan and the MS1 scan immediately
      following the current MS2 scan.
    */
    struct PuritySate_
    {
      /// Iterator pointing to the potential MS1 precursor scan
      PeakMap::ConstIterator precursorScan;
      /// Iterator pointing to the potential follow up MS1 scan
      PeakMap::ConstIterator followUpScan;

      /// Indicates if a follow up scan was found
      bool hasFollowUpScan;
      /// reference to the experiment to analyze
      const PeakMap& baseExperiment;

      /**
        @brief C'tor taking the experiment that will be analyzed.

        @param targetExp The experiment that will be analyzed.
      */
      PuritySate_(const PeakMap& targetExp);

      /**
        @brief Searches the experiment for the next MS1 spectrum with a retention time bigger then @p rt.

        @param rt The next follow up scan should have a retention bigger then this value.
      */
      void advanceFollowUp(const double rt);

      /**
        @brief Check if the currently selected follow up scan has a retention time bigger then the given value.

        @param rt The retention time to check.
      */
      bool followUpValid(const double rt) const;
    };

    /// The used quantitation method (itraq4plex, tmt6plex,..).
    const IsobaricQuantitationMethod* quant_method_;

    /// Used to select only specific types of spectra for the channel extraction.
    String selected_activation_;

    /// Allowed deviation between the expected and observed reporter ion m/z.
    Peak2D::CoordinateType reporter_mass_shift_;

    /// Minimum intensity of the precursor to be considered for quantitation.
    Peak2D::IntensityType min_precursor_intensity_;

    /// Flag if precursor with missing intensity value or missing precursor spectrum should be included or not.
    bool keep_unannotated_precursor_;

    /// Minimum reporter ion intensity to be considered for quantitation.
    Peak2D::IntensityType min_reporter_intensity_;

    /// Flag if complete quantification should be discarded if a single reporter ion has an intensity below the threshold given in IsobaricChannelExtractor::min_reporter_intensity_ .
    bool remove_low_intensity_quantifications_;

    /// Minimum precursor purity to accept the spectrum for quantitation.
    double min_precursor_purity_;

    /// Max. allowed deviation between theoretical and observed isotopic peaks of the precursor peak in the isolation window to be counted as part of the precursor.
    double max_precursor_isotope_deviation_;

    /// Flag if precursor purity will solely be computed based on the precursor scan (false), or interpolated between the precursor- and the following MS1 scan.
    bool interpolate_precursor_purity_;

    /// add channel information to the map after it has been filled
    void registerChannelsInOutputMap_(ConsensusMap& consensus_map);

    /**
      @brief Checks if the given precursor fulfills all constraints for extractions.

      @param precursor The precursor to test.
      @return $true$ if the precursor can be used for extraction, $false$ otherwise.
    */
    bool isValidPrecursor_(const Precursor& precursor) const;

    /**
      @brief Checks whether the given ConsensusFeature contains a channel that is below the given intensity threshold.

      @param cf The ConsensusFeature to check.
      @return $true$ if a low intensity reporter is contained, $false$ otherwise.
    */
    bool hasLowIntensityReporter_(const ConsensusFeature& cf) const;

    /**
      @brief Computes the purity of the precursor given an iterator pointing to the MS/MS spectrum and one to the precursor spectrum.

      @param ms2_spec Iterator pointing to the MS2 spectrum.
      @param precursor Iterator pointing to the precursor spectrum of ms2_spec.
      @return Fraction of the total intensity in the isolation window of the precursor spectrum that was assigned to the precursor.
    */
    double computePrecursorPurity_(const PeakMap::ConstIterator& ms2_spec, const PuritySate_& precursor) const;

    /**
      @brief Computes the purity of the precursor given an iterator pointing to the MS/MS spectrum and a reference to the potential precursor spectrum.

      @param ms2_spec Iterator pointing to the MS2 spectrum.
      @param precursor_spec Precursor spectrum of ms2_spec.
      @return Fraction of the total intensity in the isolation window of the precursor spectrum that was assigned to the precursor.
    */
    double computeSingleScanPrecursorPurity_(const PeakMap::ConstIterator& ms2_spec, const PeakMap::SpectrumType& precursor_spec) const;

    /**
      @brief Get the first (of potentially many) activation methods (HCD,CID,...) of this spectrum.

      @param s The spectrum
      @return Entry from Precursor::NamesOfActivationMethod or empty string.
    */
    String getActivationMethod_(const PeakMap::SpectrumType& s) const
    {
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        if (!it->getActivationMethods().empty()) return Precursor::NamesOfActivationMethod[*(it->getActivationMethods().begin())];
      }
      return "";
    }


protected:
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_() override;
  };
} // namespace

