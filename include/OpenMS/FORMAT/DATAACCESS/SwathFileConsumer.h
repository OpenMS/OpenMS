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

#ifndef OPENMS_FORMAT_DATAACCESS_SWATHFILECONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_SWATHFILECONSUMER_H

#include <boost/cast.hpp>

// Datastructures
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  /**
   * @brief Abstract base class which can consume spectra coming from SWATH experiment stored in a single file.
   *
   * The class consumes spectra which are coming from a complete SWATH experiment.
   * It expects each set of SWATH spectra to be separated by an MS1 spectrum and
   * the order of the SWATH spectra to be preserved. For example, the spectra
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
    public Interfaces::IMSDataConsumer<>
  {

public:
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    FullSwathFileConsumer() :
      ms1_counter_(0),
      ms2_counter_(0),
      ms1_map_(), // initialize to null
      consuming_possible_(true)
    {}

    ~FullSwathFileConsumer() {}

    void setExpectedSize(Size, Size) {}
    void setExperimentalSettings(const ExperimentalSettings& exp) {settings_ = exp; }

    /**
     * @brief Populate the vector of swath maps after consuming all spectra.
     *
     * Will populate the input vector with SwathMap objects which correspond to
     * the MS1 map (if present) and the MS2 maps (SWATH maps). This should be
     * called after all spectra are consumed.
     *
     * @note It is not possible to consume any more spectra after calling this
     * function (it contains finalization code and may close filestreams).
     *
     */
    void retrieveSwathMaps(std::vector<OpenSwath::SwathMap>& maps)
    {
      consuming_possible_ = false; // make consumption of further spectra / chromatograms impossble
      ensureMapsAreFilled_();
      if (ms1_map_)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(ms1_map_);
        map.lower = -1;
        map.upper = -1;
        map.ms1 = true;
        maps.push_back(map);
      }

      // TODO handle if this goes wrong ...
      assert(swath_prec_lower_.size() == swath_maps_.size());
      assert(swath_prec_upper_.size() == swath_maps_.size());

      for (Size i = 0; i < swath_maps_.size(); i++)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_maps_[i]);
        map.lower = swath_prec_lower_[i];
        map.upper = swath_prec_upper_[i];
        map.ms1 = false;
        maps.push_back(map);
      }
    }

    /// Consume a chromatogram -> should not happen when dealing with SWATH maps
    void consumeChromatogram(MapType::ChromatogramType&)
    {
      std::cerr << "Read spectrum while reading SWATH files, did not expect that!" << std::endl;
    }

    /// Consume a spectrum which may belong either to an MS1 scan or one of n MS2 (SWATH) scans
    void consumeSpectrum(MapType::SpectrumType& s)
    {
      if (!consuming_possible_)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "FullSwathFileConsumer cannot consume any more spectra after retrieveSwathMaps has been called already");
      }
      if (s.getMSLevel() == 1)
      {
        // append a new MS1 scan, set the ms2 counter to zero and proceed
        consumeMS1Spectrum_(s);
        ms2_counter_ = 0;
        ms1_counter_++;
      }
      else
      {
        // If this is the first encounter of this SWATH map, try to read the isolation windows
        if (ms2_counter_ == swath_maps_.size())
        {
          if (!s.getPrecursors().empty())
          {
            const std::vector<Precursor> prec = s.getPrecursors();
            double lower = prec[0].getIsolationWindowLowerOffset();
            double upper = prec[0].getIsolationWindowUpperOffset();
            if (prec[0].getIsolationWindowLowerOffset() > 0.0) swath_prec_lower_.push_back(lower);
            if (prec[0].getIsolationWindowUpperOffset() > 0.0) swath_prec_upper_.push_back(upper);
            swath_prec_center_.push_back(prec[0].getMZ());
          }
        }
        else if (ms2_counter_ > swath_prec_center_.size() && ms2_counter_ > swath_prec_lower_.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                           "FullSwathFileConsumer: MS2 counter is larger than size of swath maps! Are the swath_maps representing the number of read in maps?");
        }
        consumeSwathSpectrum_(s, ms2_counter_);
        ms2_counter_++;
      }
    }

protected:
    /**
     * @brief Consume an MS2 spectrum belonging to SWATH "swath_nr"
     *
     * This function should handle a Spectrum belonging to a specific SWATH
     * (indicated by swath_nr).
     *
     * @note after this call, swath_maps_.size() _must_ increase by one if
     * ms2_counter_ == swath_maps_.size() (i.e. if a new swath was encountered
     * the first time)
     */
    virtual void consumeSwathSpectrum_(MapType::SpectrumType& s, int swath_nr) = 0;
    /// @brief Consume an MS1 spectrum
    virtual void consumeMS1Spectrum_(MapType::SpectrumType& s) = 0;
    /**
     * @brief Callback function after the reading is complete
     *
     * Has to ensure that swath_maps_ and ms1_map_ are correctly populated.
     */
    virtual void ensureMapsAreFilled_() = 0;

    size_t ms1_counter_;
    size_t ms2_counter_;

    /// A list of SWATH maps and the MS1 map
    std::vector<boost::shared_ptr<MSExperiment<> > > swath_maps_;
    boost::shared_ptr<MSExperiment<> > ms1_map_;

    /// Values of lower limit, center and upper limit of the isolation windows
    std::vector<double> swath_prec_center_;
    std::vector<double> swath_prec_lower_;
    std::vector<double> swath_prec_upper_;
    /// The Experimental settings
    // (MSExperiment has no constructor using ExperimentalSettings)
    MSExperiment<> settings_;

    bool consuming_possible_;

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
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

protected:
    void addNewSwathMap_()
    {
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      swath_maps_.push_back(exp);
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, int swath_nr)
    {
      if (swath_nr == (int)swath_maps_.size())
        addNewSwathMap_();
      swath_maps_[swath_nr]->addSpectrum(s);
    }

    void addMS1Map_()
    {
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      ms1_map_ = exp;
    }

    void consumeMS1Spectrum_(MapType::SpectrumType& s)
    {
      if (!ms1_map_)
        addMS1Map_();
      ms1_map_->addSpectrum(s);
    }

    void ensureMapsAreFilled_() {}
  };

  /**
   * @brief On-disked cached implementation of FullSwathFileConsumer
   *
   * Writes all spectra immediately to disk in a user-specified caching
   * location using the CachedMzMLConsumer. Internally, it handles
   * n+1 (n SWATH + 1 MS1 map) objects of CachedMzMLConsumers which can consume the
   * spectra and write them to disk immediately.
   *
   */
  class OPENMS_DLLAPI CachedSwathFileConsumer :
    public FullSwathFileConsumer
  {

public:
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    CachedSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
      ms1_consumer_(NULL),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    ~CachedSwathFileConsumer()
    {
      // Properly delete the CachedMzMLConsumers -> free memory and _close_ filestream
      while (!swath_consumers_.empty()) {delete swath_consumers_.back(); swath_consumers_.pop_back(); }
      if (ms1_consumer_ != NULL) { delete ms1_consumer_; ms1_consumer_ = NULL; }
    }

protected:
    void addNewSwathMap_()
    {
      String meta_file = cachedir_ + basename_ + "_" + String(swath_consumers_.size()) +  ".mzML";
      String cached_file = meta_file + ".cached";
      CachedMzMLConsumer* consumer = new CachedMzMLConsumer(cached_file, true);
      consumer->setExpectedSize(nr_ms2_spectra_[swath_consumers_.size()], 0);
      swath_consumers_.push_back(consumer);

      // maps for meta data
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      swath_maps_.push_back(exp);
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, int swath_nr)
    {
      if (swath_nr == (int)swath_consumers_.size())
        addNewSwathMap_();

      swath_consumers_[swath_nr]->consumeSpectrum(s);
      swath_maps_[swath_nr]->addSpectrum(s); // append for the metadata (actual data is deleted)
    }

    void addMS1Map_()
    {
      String meta_file = cachedir_ + basename_ + "_ms1.mzML";
      String cached_file = meta_file + ".cached";
      ms1_consumer_ = new CachedMzMLConsumer(cached_file, true);
      ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      ms1_map_ = exp;
    }

    void consumeMS1Spectrum_(MapType::SpectrumType& s)
    {
      if (ms1_consumer_ == NULL)
        addMS1Map_();
      ms1_consumer_->consumeSpectrum(s);
      ms1_map_->addSpectrum(s); // append for the metadata (actual data is deleted)
    }

    void ensureMapsAreFilled_()
    {
      size_t swath_consumers_size = swath_consumers_.size();
      bool have_ms1 = (ms1_consumer_ != NULL);

      // Properly delete the CachedMzMLConsumers -> free memory and _close_ filestream
      // The filestreams to the cached data on disc can and should be closed
      // here safely. Since ensureMapsAreFilled_ is called after consuming all
      // the spectra, there will be no more spectra to append but the client
      // might already want to read after this call, so all data needs to be
      // present on disc and the filestreams closed.
      while (!swath_consumers_.empty()) { delete swath_consumers_.back(); swath_consumers_.pop_back(); }
      if (ms1_consumer_ != NULL) { delete ms1_consumer_; ms1_consumer_ = NULL; }

      if (have_ms1)
      {
        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        String meta_file = cachedir_ + basename_ + "_ms1.mzML";
        // write metadata to disk and store the correct data processing tag
        CachedmzML().writeMetadata(*ms1_map_, meta_file, true);
        MzMLFile().load(meta_file, *exp.get());
        ms1_map_ = exp;
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_consumers_size); i++)
      {
        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        String meta_file = cachedir_ + basename_ + "_" + String(i) +  ".mzML";
        // write metadata to disk and store the correct data processing tag
        CachedmzML().writeMetadata(*swath_maps_[i], meta_file, true);
        MzMLFile().load(meta_file, *exp.get());
        swath_maps_[i] = exp;
      }
    }

    CachedMzMLConsumer* ms1_consumer_;
    std::vector<CachedMzMLConsumer*> swath_consumers_;

    String cachedir_;
    String basename_;
    int nr_ms1_spectra_;
    std::vector<int> nr_ms2_spectra_;
  };
}

#endif
