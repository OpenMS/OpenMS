// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// Datastructures
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>


// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

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
    typedef PeakMap MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    FullSwathFileConsumer() :
      ms1_map_(), // initialize to null
      consuming_possible_(true),
      use_external_boundaries_(false),
      correct_window_counter_(0)
    {
      use_external_boundaries_ = !swath_map_boundaries_.empty();
    }

    /**
     * @brief Constructor
     *
     * @param swath_boundaries A vector of SwathMaps of which only the center,
     * lower and upper attributes will be used to infer the expected Swath maps.
     *
     */
    FullSwathFileConsumer(std::vector<OpenSwath::SwathMap> swath_boundaries) :
      swath_map_boundaries_(swath_boundaries),
      ms1_map_(), // initialize to null
      consuming_possible_(true),
      use_external_boundaries_(false),
      correct_window_counter_(0)
    {
      use_external_boundaries_ = !swath_map_boundaries_.empty();
    }

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
    void retrieveSwathMaps(std::vector<OpenSwath::SwathMap>& maps)
    {
      consuming_possible_ = false; // make consumption of further spectra / chromatograms impossible
      ensureMapsAreFilled_();
      if (ms1_map_)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(ms1_map_);
        map.lower = -1;
        map.upper = -1;
        map.center = -1;
        map.imLower = -1;
        map.imUpper = -1;
        map.ms1 = true;
        maps.push_back(map);
      }

      // Print warning if the lower/upper window could not be determined and we
      // required manual determination of the boundaries.
      if (!use_external_boundaries_ && correct_window_counter_ != swath_maps_.size())
      {
        std::cout << "WARNING: Could not correctly read the upper/lower limits of the SWATH windows from your input file. Read " <<
          correct_window_counter_ << " correct (non-zero) window limits (expected " << swath_maps_.size() << " windows)." << std::endl;
      }

      size_t nonempty_maps = 0;
      for (Size i = 0; i < swath_maps_.size(); i++)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_maps_[i]);
        map.lower = swath_map_boundaries_[i].lower;
        map.upper = swath_map_boundaries_[i].upper;
        map.center = swath_map_boundaries_[i].center;
        map.imLower = swath_map_boundaries_[i].imLower;
        map.imUpper = swath_map_boundaries_[i].imUpper;
        map.ms1 = false;
        maps.push_back(map);
        if (map.sptr->getNrSpectra() > 0) {nonempty_maps++;}
      }

      if (nonempty_maps != swath_map_boundaries_.size())
      {
        std::cout << "WARNING: The number nonempty maps found in the input file (" << nonempty_maps << ") is not equal to the number of provided swath window boundaries (" <<
            swath_map_boundaries_.size() << "). Please check your input." << std::endl;
      }

    }

    /// Consume a chromatogram -> should not happen when dealing with SWATH maps
    void consumeChromatogram(MapType::ChromatogramType&) override
    {
      std::cerr << "Read chromatogram while reading SWATH files, did not expect that!" << std::endl;
    }

    /**
     * @brief * Consume a spectrum which may belong either to an MS1 scan or
     * one of n MS2 (SWATH) scans
     *
     */
    void consumeSpectrum(MapType::SpectrumType& s) override
    {

      if (!consuming_possible_)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "FullSwathFileConsumer cannot consume any more spectra after retrieveSwathMaps has been called already");
      }

      if (s.getMSLevel() == 1)
      {
        consumeMS1Spectrum_(s);
      }
      else
      {
        if (s.getPrecursors().empty())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Swath scan does not provide a precursor.");
        }

        const std::vector<Precursor> prec = s.getPrecursors();
        double center = prec[0].getMZ();
        double lower = prec[0].getMZ() - prec[0].getIsolationWindowLowerOffset();
        double upper = prec[0].getMZ() + prec[0].getIsolationWindowUpperOffset();

        double lowerIm = -1; // these initial values assume IM is not present
        double upperIm = -1;

        // add IM if present
        if (s.metaValueExists("ion mobility lower limit"))
        {
          lowerIm = s.getMetaValue("ion mobility lower limit"); // want this to be -1  if no ion mobility
          upperIm = s.getMetaValue("ion mobility upper limit");
        }

        bool found = false;

        // Check if enough information is present to infer the swath
        if (center <= 0.0)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Swath scan does not provide any precursor isolation information.");
        }

        // try to match the current scan to one of the already known windows
        for (Size i = 0; i < swath_map_boundaries_.size(); i++)
        {
          // We group by the precursor mz (center of the window) since this
          // should be present in all SWATH scans.
          // also specify ion mobility, if ion mobility not present will just be -1
          if ( (std::fabs(center - swath_map_boundaries_[i].center) < 1e-6) && (std::fabs(lowerIm - swath_map_boundaries_[i].imLower) < 1e-6) && ( std::fabs(upperIm - swath_map_boundaries_[i].imUpper) < 1e-6))
          {
            found = true;
            consumeSwathSpectrum_(s, i);
            break;
          }
        }
        if (!found)
        {
          if (use_external_boundaries_)
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              String("Encountered SWATH scan with boundary ") + center + " m/z which was not present in the provided windows.");
          }
          else
          {
            consumeSwathSpectrum_(s, swath_map_boundaries_.size());

            // we found a new SWATH window
            if (lower > 0.0 && upper > 0.0)
            {correct_window_counter_++;}

            OpenSwath::SwathMap boundary;
            boundary.lower = lower;
            boundary.upper = upper;
            boundary.center = center;
            boundary.imLower = lowerIm;
            boundary.imUpper = upperIm;
            swath_map_boundaries_.push_back(boundary);

            OPENMS_LOG_DEBUG << "Adding Swath centered at " << center
              << " m/z with an isolation window of " << lower << " to " << upper
              << " m/z and IM lower limit of " <<  lowerIm << " and upper limit of " << upperIm << std::endl;
          }
        }
      }
    }

protected:

    /**
     * @brief Consume an MS2 spectrum belonging to SWATH "swath_nr"
     *
     * This function should handle a spectrum belonging to a specific SWATH
     * (indicated by swath_nr).
     *
     */
    virtual void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t swath_nr) = 0;

    /**
     * @brief Consume an MS1 spectrum
     *
     * This function should handle an MS1 spectrum.
     *
     */
    virtual void consumeMS1Spectrum_(MapType::SpectrumType& s) = 0;

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
    typedef PeakMap MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    RegularSwathFileConsumer() {}

    RegularSwathFileConsumer(std::vector<OpenSwath::SwathMap> known_window_boundaries) :
      FullSwathFileConsumer(known_window_boundaries) {}

protected:
    void addNewSwathMap_()
    {
      boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
      swath_maps_.push_back(exp);
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t swath_nr) override
    {
      while (swath_maps_.size() <= swath_nr)
      {
        addNewSwathMap_();
      }

      swath_maps_[swath_nr]->addSpectrum(s);
    }

    void addMS1Map_()
    {
      boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
      ms1_map_ = exp;
    }

    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      if (!ms1_map_)
      {
        addMS1Map_();
      }
      ms1_map_->addSpectrum(s);
    }

    void ensureMapsAreFilled_() override {}
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
    typedef PeakMap MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    CachedSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
      ms1_consumer_(nullptr),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    CachedSwathFileConsumer(std::vector<OpenSwath::SwathMap> known_window_boundaries,
            String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
      FullSwathFileConsumer(known_window_boundaries),
      ms1_consumer_(nullptr),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    ~CachedSwathFileConsumer() override
    {
      // Properly delete the MSDataCachedConsumer -> free memory and _close_ file stream
      while (!swath_consumers_.empty())
      {
        delete swath_consumers_.back();
        swath_consumers_.pop_back();
      }
      if (ms1_consumer_ != nullptr)
      {
        delete ms1_consumer_;
        ms1_consumer_ = nullptr;
      }
    }

protected:
    void addNewSwathMap_()
    {
      String meta_file = cachedir_ + basename_ + "_" + String(swath_consumers_.size()) +  ".mzML";
      String cached_file = meta_file + ".cached";
      MSDataCachedConsumer* consumer = new MSDataCachedConsumer(cached_file, true);
      consumer->setExpectedSize(nr_ms2_spectra_[swath_consumers_.size()], 0);
      swath_consumers_.push_back(consumer);

      // maps for meta data
      boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
      swath_maps_.push_back(exp);
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t swath_nr) override
    {
      while (swath_maps_.size() <= swath_nr)
      {
        addNewSwathMap_();
      }
      swath_consumers_[swath_nr]->consumeSpectrum(s); // write data to cached file; clear data from spectrum s
      swath_maps_[swath_nr]->addSpectrum(s); // append for the metadata (actual data was deleted)
    }

    void addMS1Map_()
    {
      String meta_file = cachedir_ + basename_ + "_ms1.mzML";
      String cached_file = meta_file + ".cached";
      ms1_consumer_ = new MSDataCachedConsumer(cached_file, true);
      ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
      boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
      ms1_map_ = exp;
    }

    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      if (ms1_consumer_ == nullptr)
      {
        addMS1Map_();
      }
      ms1_consumer_->consumeSpectrum(s);
      ms1_map_->addSpectrum(s); // append for the metadata (actual data is deleted)
    }

    void ensureMapsAreFilled_() override
    {
      size_t swath_consumers_size = swath_consumers_.size();
      bool have_ms1 = (ms1_consumer_ != nullptr);

      // Properly delete the MSDataCachedConsumer -> free memory and _close_ file stream
      // The file streams to the cached data on disc can and should be closed
      // here safely. Since ensureMapsAreFilled_ is called after consuming all
      // the spectra, there will be no more spectra to append but the client
      // might already want to read after this call, so all data needs to be
      // present on disc and the file streams closed.
      //
      // TODO merge with destructor code into own function!
      while (!swath_consumers_.empty())
      {
        delete swath_consumers_.back();
        swath_consumers_.pop_back();
      }
      if (ms1_consumer_ != nullptr)
      {
        delete ms1_consumer_;
        ms1_consumer_ = nullptr;
      }

      if (have_ms1)
      {
        boost::shared_ptr<PeakMap > exp(new PeakMap);
        String meta_file = cachedir_ + basename_ + "_ms1.mzML";
        // write metadata to disk and store the correct data processing tag
        Internal::CachedMzMLHandler().writeMetadata(*ms1_map_, meta_file, true);
        FileHandler().loadExperiment(meta_file, *exp.get(), {FileTypes::MZML});
        ms1_map_ = exp;
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_consumers_size); i++)
      {
        boost::shared_ptr<PeakMap > exp(new PeakMap);
        String meta_file = cachedir_ + basename_ + "_" + String(i) +  ".mzML";
        // write metadata to disk and store the correct data processing tag
        Internal::CachedMzMLHandler().writeMetadata(*swath_maps_[i], meta_file, true);
        FileHandler().loadExperiment(meta_file, *exp.get(), {FileTypes::MZML});
        swath_maps_[i] = exp;
      }
    }

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
   * Warning: no swathmaps (MS1 nor MS2) will be available when calling retrieveSwathMaps()
   *          for downstream use.
   *
   */
  class OPENMS_DLLAPI MzMLSwathFileConsumer :
    public FullSwathFileConsumer
  {

public:
    typedef PeakMap MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    MzMLSwathFileConsumer(const String& cachedir, const String& basename, Size nr_ms1_spectra, const std::vector<int>& nr_ms2_spectra) :
      ms1_consumer_(nullptr),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    MzMLSwathFileConsumer(std::vector<OpenSwath::SwathMap> known_window_boundaries,
            const String& cachedir, const String& basename, Size nr_ms1_spectra, const std::vector<int>& nr_ms2_spectra) :
      FullSwathFileConsumer(known_window_boundaries),
      ms1_consumer_(nullptr),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    ~MzMLSwathFileConsumer() override
    {
      deleteSetNull_();
    }

protected:

    void deleteSetNull_()
    {
      // Properly delete the MSDataCachedConsumer -> free memory and _close_ file stream
      while (!swath_consumers_.empty())
      {
        delete swath_consumers_.back();
        swath_consumers_.pop_back();
      }
      if (ms1_consumer_ != nullptr)
      {
        delete ms1_consumer_;
        ms1_consumer_ = nullptr;
      }
    }

    void addNewSwathMap_()
    {
      String mzml_file = cachedir_ + basename_ + "_" + String(swath_consumers_.size()) +  ".mzML";
      PlainMSDataWritingConsumer* consumer = new PlainMSDataWritingConsumer(mzml_file);
      consumer->getOptions().setCompression(true);
      consumer->setExpectedSize(nr_ms2_spectra_[swath_consumers_.size()], 0);
      swath_consumers_.push_back(consumer);
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t swath_nr) override
    {
      // only use swath_consumers_ to count how many we have already added
      while (swath_consumers_.size() <= swath_nr)
      {
        addNewSwathMap_();
      }
      swath_consumers_[swath_nr]->consumeSpectrum(s);
      s.clear(false);
    }

    void addMS1Map_()
    {
      String mzml_file = cachedir_ + basename_ + "_ms1.mzML";
      ms1_consumer_ = new PlainMSDataWritingConsumer(mzml_file);
      ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
      ms1_consumer_->getOptions().setCompression(true);
    }

    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      if (ms1_consumer_ == nullptr)
      {
        addMS1Map_();
      }
      ms1_consumer_->consumeSpectrum(s);
    }

    void ensureMapsAreFilled_() override
    {
      deleteSetNull_();
    }

    PlainMSDataWritingConsumer* ms1_consumer_;
    std::vector<PlainMSDataWritingConsumer*> swath_consumers_;

    String cachedir_;
    String basename_;
    int nr_ms1_spectra_;
    std::vector<int> nr_ms2_spectra_;
  };

}

