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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <limits>

namespace OpenMS
{

  // Aliases / chromatograms
  void MSExperiment::reserveSpaceSpectra(Size s)
  {
    spectra_.reserve(s);
  }

  void MSExperiment::reserveSpaceChromatograms(Size s)
  {
    chromatograms_.reserve(s);
  }

  /// Constructor
  MSExperiment::MSExperiment() :
    RangeManagerType(),
    ExperimentalSettings(),
    ms_levels_(),
    total_size_(0)
  {}

  /// Copy constructor
  MSExperiment::MSExperiment(const MSExperiment & source) :
    RangeManagerType(source),
    ExperimentalSettings(source),
    ms_levels_(source.ms_levels_),
    total_size_(source.total_size_),
    chromatograms_(source.chromatograms_),
    spectra_(source.spectra_)
  {}

  /// Assignment operator
  MSExperiment & MSExperiment::operator=(const MSExperiment & source)
  {
    if (&source == this) return *this;

    RangeManagerType::operator=(source);
    ExperimentalSettings::operator=(source);

    ms_levels_ = source.ms_levels_;
    total_size_ = source.total_size_;
    chromatograms_ = source.chromatograms_;
    spectra_ = source.spectra_;

    //no need to copy the alloc?!
    //alloc_

    return *this;
  }

  /// Assignment operator
  MSExperiment & MSExperiment::operator=(const ExperimentalSettings & source)
  {
    ExperimentalSettings::operator=(source);
    return *this;
  }

  /// Equality operator
  bool MSExperiment::operator==(const MSExperiment & rhs) const
  {
    return ExperimentalSettings::operator==(rhs) &&
      chromatograms_ == rhs.chromatograms_ &&
      spectra_ == rhs.spectra_;
  }

  /// Equality operator
  bool MSExperiment::operator!=(const MSExperiment & rhs) const
  {
    return !(operator==(rhs));
  }


  ///@name Iterating ranges and areas
  //@{
  /// Returns an area iterator for @p area
  MSExperiment::AreaIterator MSExperiment::areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz)
  {
    OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
    OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using AreaIterator will give invalid results!")
    //std::cout << "areaBegin: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
    return AreaIterator(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
  }

  /// Returns an invalid area iterator marking the end of an area
  MSExperiment::AreaIterator MSExperiment::areaEnd()
  {
    return AreaIterator();
  }

  /// Returns a non-mutable area iterator for @p area
  MSExperiment::ConstAreaIterator MSExperiment::areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const
  {
    OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
    OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using ConstAreaIterator will give invalid results!")
    //std::cout << "areaBeginConst: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
    return ConstAreaIterator(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
  }

  /// Returns an non-mutable invalid area iterator marking the end of an area
  MSExperiment::ConstAreaIterator MSExperiment::areaEndConst() const
  {
    return ConstAreaIterator();
  }

  /**
  @brief Fast search for spectrum range begin

  Returns the first scan which has equal or higher (>=) RT than @p rt.

  @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
  */
  MSExperiment::ConstIterator MSExperiment::RTBegin(CoordinateType rt) const
  {
    SpectrumType s;
    s.setRT(rt);
    return lower_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::RTLess());
  }

  /**
  @brief Fast search for spectrum range end (returns the past-the-end iterator)

  Returns the first scan which has higher (>) RT than @p rt.

  @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
  */
  MSExperiment::ConstIterator MSExperiment::RTEnd(CoordinateType rt) const
  {
    SpectrumType s;
    s.setRT(rt);
    return upper_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::RTLess());
  }

  /**
  @brief Fast search for spectrum range begin

  @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
  */
  MSExperiment::Iterator MSExperiment::RTBegin(CoordinateType rt)
  {
    SpectrumType s;
    s.setRT(rt);
    return lower_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::RTLess());
  }

  /**
  @brief Fast search for spectrum range end (returns the past-the-end iterator)

  @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
  */
  MSExperiment::Iterator MSExperiment::RTEnd(CoordinateType rt)
  {
    SpectrumType s;
    s.setRT(rt);
    return upper_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::RTLess());
  }

  //@}

  /**
  @name Range methods

  @note The range values (min, max, etc.) are not updated automatically. Call updateRanges() to update the values!
  */
  ///@{
  // Docu in base class
  void MSExperiment::updateRanges()
  {
    updateRanges(-1);
  }

  /**
  @brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level

  @param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
  */
  void MSExperiment::updateRanges(Int ms_level)
  {
    //clear MS levels
    ms_levels_.clear();

    //reset mz/rt/int range
    this->clearRanges();
    //reset point count
    total_size_ = 0;

    //empty
    if (spectra_.empty() && chromatograms_.empty())
    {
      return;
    }

    //update
    for (Base::iterator it = spectra_.begin(); it != spectra_.end(); ++it)
    {
      if (ms_level < Int(0) || Int(it->getMSLevel()) == ms_level)
      {
        //ms levels
        if (std::find(ms_levels_.begin(), ms_levels_.end(), it->getMSLevel()) == ms_levels_.end())
        {
          ms_levels_.push_back(it->getMSLevel());
        }

        // calculate size
        total_size_ += it->size();

        //rt
        if (it->getRT() < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(it->getRT());
        if (it->getRT() > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(it->getRT());

        //do not update mz and int when the spectrum is empty
        if (it->size() == 0) continue;

        it->updateRanges();

        //mz
        if (it->getMin()[0] < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(it->getMin()[0]);
        if (it->getMax()[0] > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(it->getMax()[0]);

        //int
        if (it->getMinInt() < RangeManagerType::int_range_.minX()) RangeManagerType::int_range_.setMinX(it->getMinInt());
        if (it->getMaxInt() > RangeManagerType::int_range_.maxX()) RangeManagerType::int_range_.setMaxX(it->getMaxInt());

      }
      // for MS level = 1 we extend the range for all the MS2 precursors
      if (ms_level == 1 && it->getMSLevel() == 2)
      {
        if (!it->getPrecursors().empty())
        {
          double pc_rt = it->getRT();
          if (pc_rt < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(pc_rt);
          if (pc_rt > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(pc_rt);
          double pc_mz = it->getPrecursors()[0].getMZ();
          if (pc_mz < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(pc_mz);
          if (pc_mz > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(pc_mz);
        }

      }

    }
    std::sort(ms_levels_.begin(), ms_levels_.end());




    if (this->chromatograms_.empty())
    {
      return;
    }

    //TODO CHROM update intensity, m/z and RT according to chromatograms as well! (done????)

    for (std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
    {

      // ignore TICs and ECs (as these are usually positioned at 0 and therefor lead to a large white margin in plots if included)
      if (it->getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM ||
        it->getChromatogramType() == ChromatogramSettings::EMISSION_CHROMATOGRAM)
      {
        continue;
      }

      // update MZ
      if (it->getMZ() < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(it->getMZ());
      if (it->getMZ() > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(it->getMZ());

      // do not update RT and in if the spectrum is empty
      if (it->size() == 0) continue;

      total_size_ += it->size();

      it->updateRanges();

      // RT
      if (it->getMin()[0] < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(it->getMin()[0]);
      if (it->getMax()[0] > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(it->getMax()[0]);

      // int
      if (it->getMinInt() < RangeManagerType::int_range_.minX()) RangeManagerType::int_range_.setMinX(it->getMinInt());
      if (it->getMaxInt() > RangeManagerType::int_range_.maxX()) RangeManagerType::int_range_.setMaxX(it->getMaxInt());
    }
  }

  /// returns the minimal m/z value
  MSExperiment::CoordinateType MSExperiment::getMinMZ() const
  {
    return RangeManagerType::pos_range_.minPosition()[1];
  }

  /// returns the maximal m/z value
  MSExperiment::CoordinateType MSExperiment::getMaxMZ() const
  {
    return RangeManagerType::pos_range_.maxPosition()[1];
  }

  /// returns the minimal retention time value
  MSExperiment::CoordinateType MSExperiment::getMinRT() const
  {
    return RangeManagerType::pos_range_.minPosition()[0];
  }

  /// returns the maximal retention time value
  MSExperiment::CoordinateType MSExperiment::getMaxRT() const
  {
    return RangeManagerType::pos_range_.maxPosition()[0];
  }

  /**
  @brief Returns RT and m/z range the data lies in.

  RT is dimension 0, m/z is dimension 1
  */
  const MSExperiment::AreaType& MSExperiment::getDataRange() const
  {
    return RangeManagerType::pos_range_;
  }

  /// returns the total number of peaks
  UInt64 MSExperiment::getSize() const
  {
    return total_size_;
  }

  /// returns an array of MS levels
  const std::vector<UInt>& MSExperiment::getMSLevels() const
  {
    return ms_levels_;
  }

  ///@}

  ///@name Sorting spectra and peaks
  ///@{
  /**
  @brief Sorts the data points by retention time

  @param sort_mz if @em true, spectra are sorted by m/z position as well
  */
  void MSExperiment::sortSpectra(bool sort_mz)
  {
    std::sort(spectra_.begin(), spectra_.end(), SpectrumType::RTLess());

    if (sort_mz)
    {
      // sort each spectrum by m/z
      for (Iterator iter = spectra_.begin(); iter != spectra_.end(); ++iter)
      {
        iter->sortByPosition();
      }
    }
  }

  /**
  @brief Sorts the data points of the chromatograms by m/z

  @param sort_rt if @em true, chromatograms are sorted by rt position as well
  */
  void MSExperiment::sortChromatograms(bool sort_rt)
  {
    // sort the chromatograms according to their product m/z
    std::sort(chromatograms_.begin(), chromatograms_.end(), ChromatogramType::MZLess());

    if (sort_rt)
    {
      for (std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
      {
        it->sortByPosition();
      }
    }
  }

  /**
  @brief Checks if all spectra are sorted with respect to ascending RT

  @param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
  */
  bool MSExperiment::isSorted(bool check_mz) const
  {
    //check RT positions
    for (Size i = 1; i < spectra_.size(); ++i)
    {
      if (spectra_[i - 1].getRT() > spectra_[i].getRT()) return false;
    }
    //check spectra
    if (check_mz)
    {
      for (Size i = 0; i < spectra_.size(); ++i)
      {
        if (!spectra_[i].isSorted()) return false;
      }
    }
    // TODO CHROM
    return true;
  }

  //@}

  /// Resets all internal values
  void MSExperiment::reset()
  {
    spectra_.clear();           //remove data
    RangeManagerType::clearRanges();           //reset range manager
    ExperimentalSettings::operator=(ExperimentalSettings());           //reset meta info
  }

  /**
  @brief Clears the meta data arrays of all contained spectra (float, integer and string arrays)

  @return @em true if meta data arrays were present and removed. @em false otherwise.
  */
  bool MSExperiment::clearMetaDataArrays()
  {
    bool meta_present = false;
    for (Size i = 0; i < spectra_.size(); ++i)
    {
      if (spectra_[i].getFloatDataArrays().size() != 0 || spectra_[i].getIntegerDataArrays().size() != 0 || spectra_[i].getStringDataArrays().size() != 0)
      {
        meta_present = true;
      }
      spectra_[i].getStringDataArrays().clear();
      spectra_[i].getIntegerDataArrays().clear();
      spectra_[i].getFloatDataArrays().clear();
    }
    return meta_present;
  }

  /// returns the meta information of this experiment (const access)
  const ExperimentalSettings& MSExperiment::getExperimentalSettings() const
  {
    return *this;
  }

  /// returns the meta information of this experiment (mutable access)
  ExperimentalSettings& MSExperiment::getExperimentalSettings()
  {
    return *this;
  }

  /// get the file path to the first MS run
  void MSExperiment::getPrimaryMSRunPath(StringList& toFill) const
  {
    std::vector<SourceFile> sfs(this->getSourceFiles());
    for (std::vector<SourceFile>::const_iterator it = sfs.begin(); it != sfs.end(); ++it)
    {
      // assemble a single location string from the URI (path to file) and file name
      String path = it->getPathToFile();
      String filename = it->getNameOfFile();

      if (path.empty() || filename.empty())
      {
        LOG_WARN << "Path or file name of primary MS run is empty. "
          << "This might be the result of incomplete conversion. "
          << "Not that tracing back e.g. identification results to the original file might more difficult." << std::endl;
      }
      else
      {
        String ms_run_location = path + "/" + filename;
        toFill.push_back(ms_run_location);
      }
    }
  }

  /**
  @brief Returns the precursor spectrum of the scan pointed to by @p iterator

  If there is no precursor scan the past-the-end iterator is returned.
  */
  MSExperiment::ConstIterator MSExperiment::getPrecursorSpectrum(ConstIterator iterator) const
  {
    if (iterator == spectra_.end() || iterator == spectra_.begin())
    {
      return spectra_.end();
    }
    UInt ms_level = iterator->getMSLevel();
    do
    {
      --iterator;
      if (iterator->getMSLevel() < ms_level)
      {
        return iterator;
      }
    } while (iterator != spectra_.begin());

    return spectra_.end();
  }

  /// Swaps the content of this map with the content of @p from
  void MSExperiment::swap(MSExperiment & from)
  {
    MSExperiment tmp;

    //swap range information
    tmp.RangeManagerType::operator=(*this);
    this->RangeManagerType::operator=(from);
    from.RangeManagerType::operator=(tmp);

    //swap experimental settings
    tmp.ExperimentalSettings::operator=(*this);
    this->ExperimentalSettings::operator=(from);
    from.ExperimentalSettings::operator=(tmp);

    // swap chromatograms
    std::swap(chromatograms_, from.chromatograms_);

    //swap peaks
    spectra_.swap(from.getSpectra());

    //swap remaining members
    ms_levels_.swap(from.ms_levels_);
    std::swap(total_size_, from.total_size_);
  }

  /// sets the spectrum list
  void MSExperiment::setSpectra(const std::vector<MSSpectrum> & spectra)
  {
    spectra_ = spectra;
  }

  /// adds a spectrum to the list
  void MSExperiment::addSpectrum(const MSSpectrum & spectrum)
  {
    spectra_.push_back(spectrum);
  }

  /// returns the spectrum list
  const std::vector<MSSpectrum>& MSExperiment::getSpectra() const
  {
    return spectra_;
  }

  /// returns the spectrum list (mutable)
  std::vector<MSSpectrum>& MSExperiment::getSpectra()
  {
    return spectra_;
  }

  /// sets the chromatogram list
  void MSExperiment::setChromatograms(const std::vector<MSChromatogram > & chromatograms)
  {
    chromatograms_ = chromatograms;
  }

  /// adds a chromatogram to the list
  void MSExperiment::addChromatogram(const MSChromatogram & chromatogram)
  {
    chromatograms_.push_back(chromatogram);
  }

  /// returns the chromatogram list
  const std::vector<MSChromatogram >& MSExperiment::getChromatograms() const
  {
    return chromatograms_;
  }

  /// returns the chromatogram list (mutable)
  std::vector<MSChromatogram >& MSExperiment::getChromatograms()
  {
    return chromatograms_;
  }

  /// @name Easy Access interface
  //@{
  /// returns a single chromatogram 
  MSChromatogram & MSExperiment::getChromatogram(Size id)
  {
    return chromatograms_[id];
  }

  /// returns a single spectrum 
  MSSpectrum & MSExperiment::getSpectrum(Size id)
  {
    return spectra_[id];
  }

  /// get the total number of spectra available
  Size MSExperiment::getNrSpectra() const
  {
    return spectra_.size();
  }

  /// get the total number of chromatograms available
  Size MSExperiment::getNrChromatograms() const
  {
    return chromatograms_.size();
  }
  //@}

  /// returns the total ion chromatogram (TIC)
  const MSChromatogram MSExperiment::getTIC() const
  {
    // The TIC is (re)calculated from the MS1 spectra. Even if MSExperiment does not contain a TIC chromatogram explicitly, it can be reported.
    MSChromatogram TIC;
    for (Base::const_iterator spec_it = spectra_.begin(); spec_it != spectra_.end(); ++spec_it)
    {
      if (spec_it->getMSLevel() == 1)
      {
        double totalIntensity = 0;
        // sum intensities of a spectrum
        for (SpectrumType::const_iterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
        {
          totalIntensity += static_cast<double>(peak_it->getIntensity());
        }
        // fill chromatogram
        ChromatogramPeakType peak;
        peak.setRT(spec_it->getRT());
        peak.setIntensity(totalIntensity);
        TIC.push_back(peak);
      }
    }
    return TIC;
  }

  /**
  @brief Clears all data and meta data

  @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
  */
  void MSExperiment::clear(bool clear_meta_data)
  {
    spectra_.clear();

    if (clear_meta_data)
    {
      clearRanges();
      this->ExperimentalSettings::operator=(ExperimentalSettings());             // no "clear" method
      chromatograms_.clear();
      ms_levels_.clear();
      total_size_ = 0;
    }
  }

  MSExperiment::SpectrumType* MSExperiment::createSpec_(PeakType::CoordinateType rt)
  {
    spectra_.emplace_back(SpectrumType());
    SpectrumType* spectrum = &(spectra_.back());
    spectrum->setRT(rt);
    spectrum->setMSLevel(1);
    return spectrum;
  }

  /*
  @brief Append a spectrum including floatdata arrays to current MSExperiment

  @param rt RT of new spectrum
  @param metadata_names Names of floatdata arrays attached to this spectrum
  @return Pointer to newly created spectrum
  */
  MSExperiment::SpectrumType* MSExperiment::createSpec_(PeakType::CoordinateType rt, const StringList& metadata_names)
  {
    SpectrumType* spectrum = createSpec_(rt);
    // create metadata arrays
    spectrum->getFloatDataArrays().reserve(metadata_names.size());
    StringList::const_iterator itm = metadata_names.begin();
    for (; itm != metadata_names.end(); ++itm)
    {
      spectrum->getFloatDataArrays().push_back(MSSpectrum::FloatDataArray());
      spectrum->getFloatDataArrays().back().setName(*itm);
    }
    return spectrum;
  }

  /// Print the contents to a stream.
  std::ostream& operator<<(std::ostream & os, const MSExperiment & exp)
  {
    os << "-- MSEXPERIMENT BEGIN --" << std::endl;

    //experimental settings
    os << static_cast<const ExperimentalSettings &>(exp);

    //spectra
    for (std::vector<MSSpectrum>::const_iterator it = exp.getSpectra().begin(); it != exp.getSpectra().end(); ++it)
    {
      os << *it;
    }

    //chromatograms
    for (std::vector<MSChromatogram >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
    {
      os << *it;
    }

    os << "-- MSEXPERIMENT END --" << std::endl;

    return os;
  }
} //namespace OpenMS
