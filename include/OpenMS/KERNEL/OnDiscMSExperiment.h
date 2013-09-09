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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_ONDISCMSEXPERIMENT_H
#define OPENMS_KERNEL_ONDISCMSEXPERIMENT_H

#include <OpenMS/INTERFACES/DataStructures.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/FORMAT/IndexedMzMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <vector>
#include <algorithm>
#include <limits>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
    @brief Representation of a mass spectrometry experiment on disk.

    @ingroup Kernel
  */
  template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
  class OPENMS_DLLAPI OnDiscMSExperiment
  {
public:

    /**
      @brief Constructor

      This initializes the object and attempts to read the indexed mzML by
      parsing the index and then reading the meta information into memory.
    */
    OnDiscMSExperiment(String filename) :
      filename_(filename),
      indexed_mzml_file_(filename)
    {
      meta_ms_experiment_ = boost::shared_ptr< MSExperiment<> >(new MSExperiment<>);

      MzMLFile f;
      PeakFileOptions options = f.getOptions();
      options.setFillData(false);
      f.setOptions(options);
      f.load(filename, *meta_ms_experiment_.get());
    }

    /// Copy constructor
    OnDiscMSExperiment(const OnDiscMSExperiment & source) :
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
    bool operator==(const OnDiscMSExperiment & rhs) const
    {
      // check if file and meta information is the same
      return filename_ == rhs.filename_ &&
        (*meta_ms_experiment_) == (*rhs.meta_ms_experiment_);
        // do not check if indexed_mzml_file_ is equal -> they have the same filename...
    }

    /// Inequality operator
    bool operator!=(const OnDiscMSExperiment & rhs) const
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

    /// alias for getSpectrum
    inline MSSpectrum<PeakT> operator[] (Size n)
    {
      return getSpectrum(n);
    }

    /**
      @brief returns a single spectrum 

      TODO: make this more efficient by reducing the copying   
    */
    MSSpectrum<PeakT> getSpectrum(Size id)
    {
      OpenMS::Interfaces::SpectrumPtr sptr = indexed_mzml_file_.getSpectrumById(id);
      MSSpectrum<PeakT> spectrum(meta_ms_experiment_->operator[](id));

      // recreate a spectrum from the data arrays!
      OpenMS::Interfaces::BinaryDataArrayPtr mz_arr = sptr->getMZArray();
      OpenMS::Interfaces::BinaryDataArrayPtr int_arr = sptr->getIntensityArray();
      spectrum.reserve(mz_arr->data.size());
      for (Size i = 0; i < mz_arr->data.size(); i++)
      {
        PeakT p;
        p.setMZ(mz_arr->data[i]);
        p.setIntensity(int_arr->data[i]);
        spectrum.push_back(p);
      }
      return spectrum;
    }

    /**
      @brief returns a single spectrum 
    */
    OpenMS::Interfaces::SpectrumPtr getSpectrumById(Size id)
    {
      return indexed_mzml_file_.getSpectrumById(id);
    }

    /**
      @brief returns a single chromatogram 

      TODO: make this more efficient by reducing the copying   
    */
    MSChromatogram<ChromatogramPeakT> getChromatogram(Size id)
    {
      OpenMS::Interfaces::ChromatogramPtr cptr = indexed_mzml_file_.getChromatogramById(id);
      MSChromatogram<ChromatogramPeakT> chromatogram(meta_ms_experiment_->getChromatogram(id));

      // recreate a chromatogram from the data arrays!
      OpenMS::Interfaces::BinaryDataArrayPtr rt_arr = cptr->getTimeArray();
      OpenMS::Interfaces::BinaryDataArrayPtr int_arr = cptr->getIntensityArray();
      chromatogram.reserve(rt_arr->data.size());
      for (Size i = 0; i < rt_arr->data.size(); i++)
      {
        ChromatogramPeakT p;
        p.setRT(rt_arr->data[i]);
        p.setIntensity(int_arr->data[i]);
        chromatogram.push_back(p);
      }

      return chromatogram;
    }

    /**
      @brief returns a single chromatogram 
    */
    OpenMS::Interfaces::ChromatogramPtr getChromatogramById(Size id)
    {
      return indexed_mzml_file_.getChromatogramById(id);
    }

private:
    /// Private Assignment operator -> we cannot copy file streams in IndexedMzMLFile
    OnDiscMSExperiment & operator=(const OnDiscMSExperiment & source) {;}

protected:

    /// The filename of the underlying data file
    const String filename_;
    /// The index of the underlying data file
    IndexedMzMLFile indexed_mzml_file_;
    /// The meta-data 
    boost::shared_ptr< MSExperiment<> > meta_ms_experiment_;
  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_ONDISCMSEXPERIMENT_H

