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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ISOBARICCHANNELEXTRACTOR_H
#define OPENMS_ANALYSIS_QUANTITATION_ISOBARICCHANNELEXTRACTOR_H

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  /**
    @brief Extracts individual channels from MS/MS spectra for isobaric labeling experiments.

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
    void extractChannels(const MSExperiment<Peak1D>& ms_exp_data, ConsensusMap& consensus_map);

private:
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
    DoubleReal min_precursor_purity_;

    /// Max. allowed deviation between theoretical and observed isotopic peaks of the precursor peak in the isolation window to be counted as part of the precursor.
    DoubleReal max_precursor_isotope_deviation_;

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

      @param ms2_spec Iterator pointing to the ms2 spectrum.
      @param precursor Iterator pointing to the precursor spectrum of ms2_spec.
      @return Fraction of the total intensity in the isolation window of the precursor spectrum that was assigned to the precursor.
    */
    DoubleReal computePrecursorPurity_(const MSExperiment<>::ConstIterator& ms2_spec, const MSExperiment<>::ConstIterator& precursor) const;

    /**
      @brief Computes the sum of all isotopic peak intensities in the window defined by (lower|upper)_mz_bound beginning from theoretical_isotope_mz.

      @param precursor Iterator pointing to the precursor spectrum used for extracting the peaks.
      @param lower_mz_bound Lower bound of the isolation window to analyze.
      @param upper_mz_bound Upper bound of the isolation window to analyze.
      @param theoretical_mz The start position for the search. Note that the intensity at this position will not included in the sum.
      @param isotope_offset The offset with which the isolation window should be searched (i.e., +/- NEUTRON_MASS/precursor_charge, +/- determines if it scans from left or right from the theoretical_isotope_mz).
    */
    DoubleReal sumPotentialIsotopePeaks_(const MSExperiment<Peak1D>::ConstIterator& precursor, const Peak1D::CoordinateType& lower_mz_bound, const Peak1D::CoordinateType& upper_mz_bound, Peak1D::CoordinateType theoretical_mz, const Peak1D::CoordinateType isotope_offset) const;

protected:
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_();
  };
} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ISOBARICCHANNELEXTRACTOR_H
