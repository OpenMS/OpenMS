// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/INTERFACES/DataStructures.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>


#include <vector>
#include <algorithm>
#include <limits>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
    @brief Representation of a mass spectrometry experiment on disk.

    @ingroup Kernel

    @note This implementation is @a not thread-safe since it keeps internally a
    single file access pointer which it moves when accessing a specific
    data item. Please provide a separate copy to each thread, e.g. 

    @code
    #pragma omp parallel for firstprivate(ondisc_map) 
    @endcode

  */
  class OPENMS_DLLAPI OnDiscMSExperiment
  {

  typedef ChromatogramPeak ChromatogramPeakT;
  typedef Peak1D PeakT;

public:

    /**
      @brief Constructor

      This initializes the object, use openFile to open a file.
    */
    OnDiscMSExperiment() {}

    /**
      @brief Open a specific file on disk.

      This tries to read the indexed mzML by parsing the index and then reading
      the meta information into memory.

      @return Whether the parsing of the file was successful (if false, the
      file most likely was not an indexed mzML file)
    */
    bool openFile(const String& filename, bool skipMetaData = false)
    {
      filename_ = filename;
      indexed_mzml_file_.openFile(filename);
      if (filename != "" && !skipMetaData)
      {
        loadMetaData_(filename);
      }
      return indexed_mzml_file_.getParsingSuccess();
    }

    /// Copy constructor
    OnDiscMSExperiment(const OnDiscMSExperiment& source) :
      filename_(source.filename_),
      indexed_mzml_file_(source.indexed_mzml_file_),
      meta_ms_experiment_(source.meta_ms_experiment_)
    {
    }

    /**
      @brief Equality operator

      This only checks whether the underlying file is the same and the parsed
      meta-information is the same. Note that the file reader (e.g. the
      std::ifstream of the file) might be in a different state.
    */
    bool operator==(const OnDiscMSExperiment& rhs) const
    {
      if (meta_ms_experiment_ == nullptr || rhs.meta_ms_experiment_ == nullptr) 
      {
        return filename_ == rhs.filename_ &&
               meta_ms_experiment_ == rhs.meta_ms_experiment_;
      }

      // check if file and meta information is the same
      return filename_ == rhs.filename_ &&
             (*meta_ms_experiment_) == (*rhs.meta_ms_experiment_);
      // do not check if indexed_mzml_file_ is equal -> they have the same filename...
    }

    /// Inequality operator
    bool operator!=(const OnDiscMSExperiment& rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @brief Checks if all spectra are sorted with respect to ascending RT

      Note that we cannot check whether all spectra are sorted (except if we
      were to load them all and check).
    */
    bool isSortedByRT() const
    {
      if (!meta_ms_experiment_) return false;

      return meta_ms_experiment_->isSorted(false);
    }

    /// alias for getNrSpectra
    inline Size size() const
    {
      return getNrSpectra();
    }

    /// returns whether spectra are empty
    inline bool empty() const
    {
      return indexed_mzml_file_.getNrSpectra() == 0;
    }

    /// get the total number of spectra available
    inline Size getNrSpectra() const
    {
      return indexed_mzml_file_.getNrSpectra();
    }

    /// get the total number of chromatograms available
    inline Size getNrChromatograms() const
    {
      return indexed_mzml_file_.getNrChromatograms();
    }

    /// returns the meta information of this experiment (const access)
    boost::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const
    {
      return boost::static_pointer_cast<const ExperimentalSettings>(meta_ms_experiment_);
    }

    boost::shared_ptr<PeakMap> getMetaData() const
    {
      return meta_ms_experiment_;
    }

    /// alias for getSpectrum
    inline MSSpectrum operator[](Size n)
    {
      return getSpectrum(n);
    }

    /**
      @brief returns a single spectrum

      @param id The index of the spectrum
    */
    MSSpectrum getSpectrum(Size id)
    {
      if (!meta_ms_experiment_) return indexed_mzml_file_.getMSSpectrumById(int(id));

      MSSpectrum spectrum(meta_ms_experiment_->operator[](id));
      indexed_mzml_file_.getMSSpectrumById(int(id), spectrum);
      return spectrum;
    }

    /**
      @brief returns a single spectrum
    */
    OpenMS::Interfaces::SpectrumPtr getSpectrumById(Size id)
    {
      return indexed_mzml_file_.getSpectrumById((int)id);
    }

    /**
      @brief returns a single chromatogram

      @param id The index of the chromatogram
    */
    MSChromatogram getChromatogram(Size id)
    {
      if (!meta_ms_experiment_) return indexed_mzml_file_.getMSChromatogramById(int(id));

      MSChromatogram chromatogram(meta_ms_experiment_->getChromatogram(id));
      indexed_mzml_file_.getMSChromatogramById(int(id), chromatogram);
      return chromatogram;
    }

    /**
      @brief returns a single chromatogram

      @param id The native identifier of the chromatogram
    */
    MSChromatogram getChromatogramByNativeId(const std::string& id);

    /**
      @brief returns a single spectrum

      @param id The native identifier of the spectrum
    */
    MSSpectrum getSpectrumByNativeId(const std::string& id);

    /**
      @brief returns a single chromatogram
    */
    OpenMS::Interfaces::ChromatogramPtr getChromatogramById(Size id)
    {
      return indexed_mzml_file_.getChromatogramById(id);
    }

    /// sets whether to skip some XML checks and be fast instead
    void setSkipXMLChecks(bool skip)
    {
      indexed_mzml_file_.setSkipXMLChecks(skip);
    }

private:

    /// Private Assignment operator -> we cannot copy file streams in IndexedMzMLHandler
    OnDiscMSExperiment& operator=(const OnDiscMSExperiment& /* source */);

    void loadMetaData_(const String& filename);

    MSChromatogram getMetaChromatogramById_(const std::string& id);

    MSSpectrum getMetaSpectrumById_(const std::string& id);

protected:

    /// The filename of the underlying data file
    String filename_;
    /// The index of the underlying data file
    Internal::IndexedMzMLHandler indexed_mzml_file_;
    /// The meta-data
    boost::shared_ptr<PeakMap> meta_ms_experiment_;
    /// Mapping of chromatogram native ids to offsets
    std::unordered_map< std::string, Size > chromatograms_native_ids_;
    /// Mapping of spectra native ids to offsets
    std::unordered_map< std::string, Size > spectra_native_ids_;
  };

typedef OpenMS::OnDiscMSExperiment OnDiscPeakMap;

} // namespace OpenMS


