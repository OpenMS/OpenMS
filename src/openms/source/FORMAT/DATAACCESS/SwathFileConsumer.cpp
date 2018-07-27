// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>

namespace OpenMS
{

  FullSwathFileConsumer::FullSwathFileConsumer() :
    ms1_map_(), // initialize to null
    consuming_possible_(true),
    use_external_boundaries_(false),
    correct_window_counter_(0)
  {
    use_external_boundaries_ = !swath_map_boundaries_.empty();
  }

  FullSwathFileConsumer::FullSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& swath_boundaries) :
    swath_map_boundaries_(swath_boundaries),
    ms1_map_(), // initialize to null
    consuming_possible_(true),
    use_external_boundaries_(false),
    correct_window_counter_(0)
  {
    use_external_boundaries_ = !swath_map_boundaries_.empty();
  }

  void FullSwathFileConsumer::retrieveSwathMaps(std::vector<OpenSwath::SwathMap>& maps)
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

  void FullSwathFileConsumer::consumeChromatogram(ChromatogramType&)
  {
    std::cerr << "Read chromatogram while reading SWATH files, did not expect that!" << std::endl;
  }

  void FullSwathFileConsumer::consumeSpectrum(SpectrumType& s)
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
        if (std::fabs(center - swath_map_boundaries_[i].center) < 1e-6)
        {
          found = true;
          consumeSwathSpectrum_(s, i);
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
          swath_map_boundaries_.push_back(boundary);

          LOG_DEBUG << "Adding Swath centered at " << center
            << " m/z with an isolation window of " << lower << " to " << upper
            << " m/z." << std::endl;
        }
      }
    }
  }

  RegularSwathFileConsumer::RegularSwathFileConsumer() {}

  RegularSwathFileConsumer::RegularSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries) :
    FullSwathFileConsumer(known_window_boundaries) {}

  void RegularSwathFileConsumer::addNewSwathMap_()
  {
    boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
    swath_maps_.push_back(exp);
  }

  void RegularSwathFileConsumer::consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr) 
  {
    while (swath_maps_.size() <= swath_nr)
    {
      addNewSwathMap_();
    }

    swath_maps_[swath_nr]->addSpectrum(s);
  }

  void RegularSwathFileConsumer::addMS1Map_()
  {
    boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
    ms1_map_ = exp;
  }

  void RegularSwathFileConsumer::consumeMS1Spectrum_(SpectrumType& s) 
  {
    if (!ms1_map_)
    {
      addMS1Map_();
    }
    ms1_map_->addSpectrum(s);
  }

  void RegularSwathFileConsumer::ensureMapsAreFilled_() {}




  CachedSwathFileConsumer::CachedSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
    ms1_consumer_(nullptr),
    swath_consumers_(),
    cachedir_(cachedir),
    basename_(basename),
    nr_ms1_spectra_(nr_ms1_spectra),
    nr_ms2_spectra_(nr_ms2_spectra)
  {}

  CachedSwathFileConsumer::CachedSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries,
      String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
    FullSwathFileConsumer(known_window_boundaries),
    ms1_consumer_(nullptr),
    swath_consumers_(),
    cachedir_(cachedir),
    basename_(basename),
    nr_ms1_spectra_(nr_ms1_spectra),
    nr_ms2_spectra_(nr_ms2_spectra)
  {}

  CachedSwathFileConsumer::~CachedSwathFileConsumer()
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

  void CachedSwathFileConsumer::addNewSwathMap_()
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

  void CachedSwathFileConsumer::consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr)
  {
    while (swath_maps_.size() <= swath_nr)
    {
      addNewSwathMap_();
    }
    swath_consumers_[swath_nr]->consumeSpectrum(s);
    swath_maps_[swath_nr]->addSpectrum(s); // append for the metadata (actual data is deleted)
  }

  void CachedSwathFileConsumer::addMS1Map_()
  {
    String meta_file = cachedir_ + basename_ + "_ms1.mzML";
    String cached_file = meta_file + ".cached";
    ms1_consumer_ = new MSDataCachedConsumer(cached_file, true);
    ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
    boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
    ms1_map_ = exp;
  }

  void CachedSwathFileConsumer::consumeMS1Spectrum_(SpectrumType& s)
  {
    if (ms1_consumer_ == nullptr)
    {
      addMS1Map_();
    }
    ms1_consumer_->consumeSpectrum(s);
    ms1_map_->addSpectrum(s); // append for the metadata (actual data is deleted)
  }

  void CachedSwathFileConsumer::ensureMapsAreFilled_()
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
      MzMLFile().load(meta_file, *exp.get());
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
      MzMLFile().load(meta_file, *exp.get());
      swath_maps_[i] = exp;
    }
  }

  MzMLSwathFileConsumer::MzMLSwathFileConsumer(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
      ms1_consumer_(nullptr),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
  {}

  MzMLSwathFileConsumer::MzMLSwathFileConsumer(const std::vector<OpenSwath::SwathMap>& known_window_boundaries,
      String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
    FullSwathFileConsumer(known_window_boundaries),
    ms1_consumer_(nullptr),
    swath_consumers_(),
    cachedir_(cachedir),
    basename_(basename),
    nr_ms1_spectra_(nr_ms1_spectra),
    nr_ms2_spectra_(nr_ms2_spectra)
  {}

  MzMLSwathFileConsumer::~MzMLSwathFileConsumer()
  {
    deleteSetNull_();
  }

  void MzMLSwathFileConsumer::deleteSetNull_()
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

  void MzMLSwathFileConsumer::addNewSwathMap_()
  {
    String mzml_file = cachedir_ + basename_ + "_" + String(swath_consumers_.size()) +  ".mzML";
    PlainMSDataWritingConsumer* consumer = new PlainMSDataWritingConsumer(mzml_file);
    consumer->getOptions().setCompression(true);
    consumer->setExpectedSize(nr_ms2_spectra_[swath_consumers_.size()], 0);
    swath_consumers_.push_back(consumer);
  }

  void MzMLSwathFileConsumer::consumeSwathSpectrum_(SpectrumType& s, size_t swath_nr)
  {
    // only use swath_maps_ to count how many we have already added
    while (swath_consumers_.size() <= swath_nr)
    {
      addNewSwathMap_();
    }
    swath_consumers_[swath_nr]->consumeSpectrum(s);
    s.clear(false);
  }

  void MzMLSwathFileConsumer::addMS1Map_()
  {
    String mzml_file = cachedir_ + basename_ + "_ms1.mzML";
    ms1_consumer_ = new PlainMSDataWritingConsumer(mzml_file);
    ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
    ms1_consumer_->getOptions().setCompression(true);
    boost::shared_ptr<PeakMap > exp(new PeakMap(settings_));
    // ms1_map_ = exp;
  }

  void MzMLSwathFileConsumer::consumeMS1Spectrum_(SpectrumType& s)
  {
    if (ms1_consumer_ == nullptr)
    {
      addMS1Map_();
    }
    ms1_consumer_->consumeSpectrum(s);
    s.clear(false);
  }

  void MzMLSwathFileConsumer::ensureMapsAreFilled_()
  {
    deleteSetNull_();
  }

} // namespace OpenMS
