// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <boost/cast.hpp>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenSwath
{
  class SwathMap;
}

namespace OpenMS
{
  class MSDataCachedConsumer;
  class PlainMSDataWritingConsumer;

  /**
   * @brief Abstract base class which can consume spectra coming from SWATH experiment stored in a single file.
   *
   * The class consumes spectra which are coming from a complete SWATH
   * experiment. It will group MS2 spectra by their precursor m/z, assuming
   * that they correspond to the same SWATH window.  For example, the spectra
   * could be arranged in the following fashion:
   *
   * - MS1 Spectrum (no precursor)
   * - MS2 Spectrum (precursor = [400,425])
   * - MS2 Spectrum (precursor = [425,450])
   * - [...]
   * - MS2 Spectrum (precursor = [1175,1200])
   * - MS1 Spectrum (no precursor)
   * - MS2 Spectrum (precursor = [400,425])
   * - MS2 Spectrum (precursor = [425,450])
   * - [...]
   *
   * Base classes are expected to implement functions consuming a spectrum coming
   * from a specific SWATH or an MS1 spectrum and a final function
   * ensureMapsAreFilled_ after which the swath_maps_ vector needs to contain
   * valid pointers to MSExperiment.
   *
   * In addition it is possible to provide the swath boundaries and the read in
   * spectra will be matched by their precursor m/z to the "center" attribute
   * of the provided Swath maps.
   *
   * Usage:
   *
   * @code
   * FullSwathFileConsumer * dataConsumer;
   * // assign dataConsumer to an implementation of FullSwathFileConsumer
   * MzMLFile().transform(file, dataConsumer);
   * dataConsumer->retrieveSwathMaps(maps);
   * @endcode
   *
   */
  class OPENMS_DLLAPI FullSwathFileConsumer :
    public Interfaces::IMSDataConsumer
  {

public:
    typedef MSSpectrum SpectrumType;
    typedef MSChromatogram ChromatogramType;

    FullSwathFileConsumer();

    /**
     * @brief Constructor
     *
     * @param swath_boundaries A vector of SwathMaps of which only the center,
     * lower and upper attributes will be used to infer the expected Swath maps.
     *
     */
    FullSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& swath_boundaries);

    ~FullSwathFileConsumer() override {}

    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings& exp) override {settings_ = exp; }

    /**
     * @brief Populate the vector of swath maps after consuming all spectra.
     *
     * Will populate the input vector with SwathMap objects which correspond to
     * the MS1 map (if present) and the MS2 maps (SWATH maps). This should be
     * called after all spectra are consumed.
     *
     * @note It is not possible to consume any more spectra after calling this
     * function (it contains finalization code and may close file streams).
     *
     */
    void retrieveSwathMaps(std::vector<OpenSwath::SwathMap>& maps);

    /// Consume a chromatogram -> should not happen when dealing with SWATH maps
    void consumeChromatogram(ChromatogramType&) override;

    /**
     * @brief * Consume a spectrum which may belong either to an MS1 scan or
     * one of n MS2 (SWATH) scans
     *
     */
    void consumeSpectrum(SpectrumType& s) override;

protected:

    /**
     * @brief Consume an MS2 spectrum belonging to SWATH "swath_nr"
     *
     * This function should handle a spectrum belonging to a specific SWATH
     * (indicated by swath_nr).
     *
     */
    virtual void consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr) = 0;

    /**
     * @brief Consume an MS1 spectrum
     *
     * This function should handle an MS1 spectrum.
     *
     */
    virtual void consumeMS1Spectrum_(SpectrumType& s) = 0;

    /**
     * @brief Callback function after the reading is complete
     *
     * Has to ensure that swath_maps_ and ms1_map_ are correctly populated.
     */
    virtual void ensureMapsAreFilled_() = 0;

    /// A list of Swath map identifiers (lower/upper boundary and center)
    std::vector<OpenSwath::SwathMap> swath_map_boundaries_;

    /// A list of SWATH maps and the MS1 map
    std::vector<boost::shared_ptr<PeakMap > > swath_maps_;
    boost::shared_ptr<PeakMap > ms1_map_;

    /// The Experimental settings
    // (MSExperiment has no constructor using ExperimentalSettings)
    PeakMap settings_;

    /// Whether further spectra can still be consumed
    bool consuming_possible_;

    /// Whether to use external input for SWATH boundaries
    bool use_external_boundaries_;

    /// How many windows were correctly annotated (non-zero window limits)
    size_t correct_window_counter_;

  };

  /**
   * @brief In-memory implementation of FullSwathFileConsumer
   *
   * Keeps all the spectra in memory by just appending them to an MSExperiment.
   *
   */
  class OPENMS_DLLAPI RegularSwathFileConsumer :
    public FullSwathFileConsumer
  {

public:
    typedef MSSpectrum SpectrumType;
    typedef MSChromatogram ChromatogramType;

    RegularSwathFileConsumer();

    RegularSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries);

protected:

    void addNewSwathMap_();

    void consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr) override;

    void addMS1Map_();

    void consumeMS1Spectrum_(SpectrumType& s) override;

    void ensureMapsAreFilled_() override;
  };

  /**
   * @brief On-disk cached implementation of FullSwathFileConsumer
   *
   * Writes all spectra immediately to disk in a user-specified caching
   * location using the MSDataCachedConsumer. Internally, it handles
   * n+1 (n SWATH + 1 MS1 map) objects of MSDataCachedConsumer which can consume the
   * spectra and write them to disk immediately.
   *
   */
  class OPENMS_DLLAPI CachedSwathFileConsumer :
    public FullSwathFileConsumer
  {

public:
    typedef MSSpectrum SpectrumType;
    typedef MSChromatogram ChromatogramType;

    CachedSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra);

    CachedSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries,
            String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra);

    ~CachedSwathFileConsumer() override;

protected:
    void addNewSwathMap_();

    void consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr) override;

    void addMS1Map_();

    void consumeMS1Spectrum_(SpectrumType& s) override;

    void ensureMapsAreFilled_() override;

    MSDataCachedConsumer* ms1_consumer_;
    std::vector<MSDataCachedConsumer*> swath_consumers_;

    String cachedir_;
    String basename_;
    int nr_ms1_spectra_;
    std::vector<int> nr_ms2_spectra_;
  };

  /**
   * @brief On-disk mzML implementation of FullSwathFileConsumer
   *
   * Writes all spectra immediately to disk to an mzML file location using the
   * PlainMSDataWritingConsumer. Internally, it handles n+1 (n SWATH + 1 MS1
   * map) objects of MSDataCachedConsumer which can consume the spectra and
   * write them to disk immediately.
   *
   */
  class OPENMS_DLLAPI MzMLSwathFileConsumer :
    public FullSwathFileConsumer
  {

public:
    typedef MSSpectrum SpectrumType;
    typedef MSChromatogram ChromatogramType;

    MzMLSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra);

    MzMLSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries,
            String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra);

    ~MzMLSwathFileConsumer() override;

protected:

    void deleteSetNull_();

    void addNewSwathMap_();

    void consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr) override;

    void addMS1Map_();

    void consumeMS1Spectrum_(SpectrumType& s) override;

    void ensureMapsAreFilled_() override;

    PlainMSDataWritingConsumer* ms1_consumer_;
    std::vector<PlainMSDataWritingConsumer*> swath_consumers_;

    String cachedir_;
    String basename_;
    int nr_ms1_spectra_;
    std::vector<int> nr_ms2_spectra_;
  };

}

