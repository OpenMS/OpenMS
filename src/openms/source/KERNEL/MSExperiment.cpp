// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <limits>
#include <unordered_set>

namespace OpenMS
{
  /// Constructor
  MSExperiment::MSExperiment() :
    RangeManagerContainerType(),
    ExperimentalSettings()
  {}

  /// Copy constructor
  MSExperiment::MSExperiment(const MSExperiment & source) = default;

  /// Assignment operator
  MSExperiment & MSExperiment::operator=(const MSExperiment & source)
  {
    if (&source == this)
    {
      return *this;
    }
    RangeManagerContainerType::operator=(source);
    ExperimentalSettings::operator=(source);

    chromatograms_ = source.chromatograms_;
    spectra_ = source.spectra_;

    //no need to copy the alloc?!
    //alloc_

    return *this;
  }

  /// Assignment operator
  MSExperiment& MSExperiment::operator=(const ExperimentalSettings & source)
  {
    ExperimentalSettings::operator=(source);
    return *this;
  }

  MSExperiment::~MSExperiment() = default;

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

  void MSExperiment::reserveSpaceSpectra(Size s)
  {
    spectra_.reserve(s);
  }

  void MSExperiment::reserveSpaceChromatograms(Size s)
  {
    chromatograms_.reserve(s);
  }

  ///@name Iterating ranges and areas
  //@{
  /// Returns an area iterator for @p area
  MSExperiment::AreaIterator MSExperiment::areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt ms_level)
  {
    OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
    OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using AreaIterator will give invalid results!")
    auto [min_im, max_im] = RangeMobility{}.getNonEmptyRange(); // a full range
    auto area = AreaIterator::Param(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), ms_level);
    area.lowMZ(min_mz).highMZ(max_mz).lowIM(min_im).highIM(max_im);
    return AreaIterator(area);
  }

  MSExperiment::AreaIterator MSExperiment::areaBegin(const RangeManagerType& range, UInt ms_level)
  {
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using ConstAreaIterator will give invalid results!")
    auto [min_rt, max_rt] = range.RangeRT::getNonEmptyRange();
    auto [min_mz, max_mz] = range.RangeMZ::getNonEmptyRange();
    auto [min_im, max_im] = range.RangeMobility::getNonEmptyRange();
    auto area = AreaIterator::Param(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), ms_level);
    area.lowMZ(min_mz).highMZ(max_mz).lowIM(min_im).highIM(max_im);
    return AreaIterator(area);
  }

  /// Returns an invalid area iterator marking the end of an area
  MSExperiment::AreaIterator MSExperiment::areaEnd()
  {
    return AreaIterator();
  }

  /// Returns a non-mutable area iterator for @p area
  MSExperiment::ConstAreaIterator MSExperiment::areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt ms_level) const
  {
    OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
    OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using ConstAreaIterator will give invalid results!")
    auto [min_im, max_im] = RangeMobility{}.getNonEmptyRange(); // a full range
    auto area = ConstAreaIterator::Param(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), ms_level);
    area.lowMZ(min_mz).highMZ(max_mz).lowIM(min_im).highIM(max_im);
    return ConstAreaIterator(area);
  }

  MSExperiment::ConstAreaIterator MSExperiment::areaBeginConst(const RangeManagerType& range, UInt ms_level) const
  {
    OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using ConstAreaIterator will give invalid results!")
    auto [min_rt, max_rt] = range.RangeRT::getNonEmptyRange();
    auto [min_mz, max_mz] = range.RangeMZ::getNonEmptyRange();
    auto [min_im, max_im] = range.RangeMobility::getNonEmptyRange();
    auto area = ConstAreaIterator::Param(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), ms_level);
    area.lowMZ(min_mz).highMZ(max_mz).lowIM(min_im).highIM(max_im);
    return ConstAreaIterator(area);
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

  MSExperiment::ConstIterator MSExperiment::IMBegin(CoordinateType im) const
  {
    SpectrumType s;
    s.setDriftTime(im);
    return lower_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::IMLess());
  }

  MSExperiment::ConstIterator MSExperiment::IMEnd(CoordinateType im) const
  {
    SpectrumType s;
    s.setDriftTime(im);
    return upper_bound(spectra_.begin(), spectra_.end(), s, SpectrumType::IMLess());
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
  @brief Updates the m/z, intensity, retention time, ion mobility and MS level ranges of all spectra with a certain ms level

  @param ms_level MS level to consider for m/z range, RT range, intensity range and ion mobility (if negative, all MS levels are used)
  */
  void MSExperiment::updateRanges(Int ms_level)
  {
    // reset mz/rt/int range
    this->clearRanges();

    // empty
    if (spectra_.empty() && chromatograms_.empty())
    {
      return;
    }

    // update
    for (Base::iterator it = spectra_.begin(); it != spectra_.end(); ++it)
    {
      if (ms_level < Int(0) || Int(it->getMSLevel()) == ms_level)
      {
        // ranges
        this->extendRT(it->getRT()); // RT
        // m/z, intensity and ion mobility from spectrum's range
        it->updateRanges();
        this->extend(*it);
      }
      // for MS level = 1 we extend the range for all the MS2 precursors
      if (ms_level == 1 && it->getMSLevel() == 2)
      {
        if (!it->getPrecursors().empty())
        {
          this->extendRT(it->getRT());
          this->extendMZ(it->getPrecursors()[0].getMZ());
        }
      }
    }

    if (this->chromatograms_.empty())
    {
      return;
    }

    // update intensity, m/z and RT according to chromatograms as well:
    for (ChromatogramType& cp : chromatograms_)
    {
      // update range of EACH chrom, if we need them individually later
      cp.updateRanges();

      // ignore TICs and ECs for the whole experiments range (as these are usually positioned at 0 and therefor lead to a large white margin in plots if included)
      if (cp.getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM ||
        cp.getChromatogramType() == ChromatogramSettings::EMISSION_CHROMATOGRAM)
      {
        continue;
      }

      // ranges
      this->extendMZ(cp.getMZ());// MZ
      this->extend(cp);// RT and intensity from chroms's range
    }
  }

  /// returns the total number of peaks
  UInt64 MSExperiment::getSize() const
  {    
    Size total_size{};
    for (const auto& spec : spectra_) total_size += spec.size(); // sum up all peaks in all spectra
    for (const auto& chrom : chromatograms_) total_size += chrom.size(); // sum up all peaks in all chromatograms
    return total_size;
  }

  /// returns an array of MS levels (calculated on demand)
  std::vector<UInt> MSExperiment::getMSLevels() const
  {
    std::unordered_set<UInt> level_set;
    for (const auto& spec : spectra_)
    {
      level_set.insert(spec.getMSLevel());
    }
    
    std::vector<UInt> ms_levels(level_set.begin(), level_set.end());
    std::sort(ms_levels.begin(), ms_levels.end());
    return ms_levels;
  }

  const String sqMassRunID = "sqMassRunID";

  UInt64 MSExperiment::getSqlRunID() const
  {
    if (metaValueExists(sqMassRunID))
    {
      return getMetaValue(sqMassRunID);
    }
    return 0;
  }

  void MSExperiment::setSqlRunID(UInt64 id)
  {
    setMetaValue(sqMassRunID, id);
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
      for (ChromatogramType& cp : chromatograms_)
      {
        cp.sortByPosition();
      }
    }
  }

  /**
  @brief Checks if all spectra are sorted with respect to ascending RT

  @param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
  */
  bool MSExperiment::isSorted(bool check_mz) const
  {
    // check RT positions
    for (Size i = 1; i < spectra_.size(); ++i)
    {
      if (spectra_[i - 1].getRT() > spectra_[i].getRT())
      {
        return false;
      }
    }
    // check spectra
    if (check_mz)
    {
      for (Size i = 0; i < spectra_.size(); ++i)
      {
        if (!spectra_[i].isSorted())
        {
          return false;
        }
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
      if (!spectra_[i].getFloatDataArrays().empty() 
        || !spectra_[i].getIntegerDataArrays().empty() 
        || !spectra_[i].getStringDataArrays().empty())
      {
        meta_present = true;
      }
      spectra_[i].getStringDataArrays().clear();
      spectra_[i].getStringDataArrays().shrink_to_fit();
      spectra_[i].getIntegerDataArrays().clear();
      spectra_[i].getIntegerDataArrays().shrink_to_fit();
      spectra_[i].getFloatDataArrays().clear();
      spectra_[i].getFloatDataArrays().shrink_to_fit();
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
    for (const SourceFile& ss : sfs)
    {
      // assemble a single location string from the URI (path to file) and file name
      String path = ss.getPathToFile();
      String filename = ss.getNameOfFile();

      if (path.empty() || filename.empty())
      {
        OPENMS_LOG_WARN << "Path or file name of primary MS run is empty. "
          << "This might be the result of incomplete conversion. "
          << "Not that tracing back e.g. identification results to the original file might more difficult." << std::endl;
      }
      else
      {
        // use Windows or UNIX path separator?
        String actual_path = path.hasPrefix("file:///") ? path.substr(8) : path;
        String sep = (actual_path.has('\\') && !actual_path.has('/')) ? "\\" : "/";
        String ms_run_location = path + sep + filename;
        toFill.push_back(ms_run_location);
      }
    }
  }

  /**
  @brief Returns the precursor spectrum of the scan pointed to by @p iterator

  If there is no precursor scan the past-the-end iterator is returned.
  This assumes that precursors occur somewhere before the current spectrum
  but not necessarily the first one from the last MS level (we double-check with
  the annotated precursorList.
  */
  MSExperiment::ConstIterator MSExperiment::getPrecursorSpectrum(ConstIterator iterator) const
  {
    // if we are after the end or at the beginning where we can't go "up"
    if (iterator == spectra_.end() || iterator == spectra_.begin())
    {
      return spectra_.end();
    }
    UInt ms_level = iterator->getMSLevel();

    if (ms_level == 1) // assumes there is not level 0
    {
      return spectra_.end();
    }

    if (!iterator->getPrecursors().empty())
    {
      //TODO warn about taking first with the blocking LOG_WARN in such a central class?
      //if (iterator->getPrecursors().size() > 1) ...

      const auto precursor = iterator->getPrecursors()[0];
      if (precursor.metaValueExists("spectrum_ref"))
      {
        String ref = precursor.getMetaValue("spectrum_ref");
        auto tmp_spec_iter = iterator; // such that we can reiterate later
        do
        {
          --tmp_spec_iter;
          if ((ms_level - tmp_spec_iter->getMSLevel() == 1) && (tmp_spec_iter->getNativeID() == ref))
          {
            return tmp_spec_iter;
          }
        } while (tmp_spec_iter != spectra_.begin());
      }
    }

    // if no precursor annotation was found or it did not have a spectrum reference,
    // just
    do
    {
      --iterator;
      if (ms_level - iterator->getMSLevel() == 1)
      {
        return iterator;
      }
    } while (iterator != spectra_.begin());

    return spectra_.end();
  }

  // same as above but easier to wrap in python
  int MSExperiment::getPrecursorSpectrum(int zero_based_index) const
  {
    auto spec = spectra_.cbegin();
    spec += zero_based_index;
    auto pc_spec = getPrecursorSpectrum(spec);
    if (pc_spec == spectra_.cend()) return -1;
    return pc_spec - spectra_.cbegin(); 
  }

  MSExperiment::ConstIterator MSExperiment::getFirstProductSpectrum(ConstIterator parent_iterator) const
  {
    // if we are already at or after the end -> there can be no product after it
    if (parent_iterator == spectra_.end() 
      || parent_iterator == spectra_.end() - 1)
    {
      return spectra_.end();
    }

    UInt parent_ms_level = parent_iterator->getMSLevel();
    const auto& parent_native_id = parent_iterator->getNativeID();

    auto it = parent_iterator; // such that we can reiterate later
    it++; // start at the next spectrum

    while (it != spectra_.end())
    { 
      if (it->getMSLevel() < parent_ms_level) return spectra_.end();

      if ((it->getMSLevel() - parent_ms_level) == 1) 
      { // it is a potential product spectrum (one level higher than parent)

        // does it have precursors referencing the parents?
        if (it->getPrecursors().empty()) 
        {
          ++it; // no precursors, so we can't check if it is a product of the parent
          continue;
        }      

        // warn if there are multiple precursors (should not happen)
        if (it->getPrecursors().size() > 1)
        {
            OPENMS_LOG_WARN << "Spectrum at index " << std::distance(spectra_.begin(), it)
                      << " has multiple precursors. Only the first precursor will be considered."
                      << std::endl;
        }

        // check if it has the parent a precursor
        const auto precursor = it->getPrecursors()[0];
        String ref = precursor.getMetaValue("spectrum_ref", "");  
        if (!ref.empty() && ref == parent_native_id)
        {
          return it;
        }     
      }
      ++it; 
    } 

    return spectra_.end();
  }

  // same as above but easier to wrap in python
  int MSExperiment::getFirstProductSpectrum(int zero_based_index) const
  {
    if (zero_based_index < 0 
      || zero_based_index >= static_cast<int>(spectra_.size()))
    {
      return -1;
    }
    auto spec = spectra_.cbegin();
    spec += zero_based_index;
    auto pc_spec = getFirstProductSpectrum(spec);
    if (pc_spec == spectra_.cend()) return -1;
    return pc_spec - spectra_.cbegin();
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

  }

  /// sets the spectrum list
  void MSExperiment::setSpectra(const std::vector<MSSpectrum> & spectra)
  {
    spectra_ = spectra;
  }

  void MSExperiment::setSpectra(std::vector<MSSpectrum> && spectra)
  {
    spectra_ = std::move(spectra);
  }

  /// adds a spectrum to the list
  void MSExperiment::addSpectrum(const MSSpectrum & spectrum)
  {
    spectra_.push_back(spectrum);
  }

  void MSExperiment::addSpectrum(MSSpectrum && spectrum)
  {
    spectra_.push_back(std::move(spectrum));
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

  /// Returns the closest(=nearest) spectrum in retention time to the given RT
  MSExperiment::ConstIterator MSExperiment::getClosestSpectrumInRT(const double RT) const
  {
    auto above = RTBegin(RT);           // the spec above or equal to our RT
    if (above == begin()) return above; // we hit the first element, or no spectra (begin==end)
    if (above == end()) return --above; // queried beyond last spec, but we know there are spectra, so `--above` is safe
    // we are between two spectra
    auto diff_left = RT - (above - 1)->getRT();
    auto diff_right = above->getRT() - RT;
    if (diff_left < diff_right) --above;
    return above;
  }
  MSExperiment::Iterator MSExperiment::getClosestSpectrumInRT(const double RT)
  {
    return begin() + std::distance(cbegin(), const_cast<const MSExperiment*>(this)->getClosestSpectrumInRT(RT));
  }

  /// Returns the closest(=nearest) spectrum in retention time to the given RT of a certain MS level
  MSExperiment::ConstIterator MSExperiment::getClosestSpectrumInRT(const double RT, UInt ms_level) const
  {
    auto above = RTBegin(RT); // the spec above or equal to our RT
    auto below = above; // for later
    // search for the next available spec to the right with correct MS level
    while (above != end() && above->getMSLevel() != ms_level)
    {
      ++above;
    }
    if (above == begin()) return above; // we hit the first element; or no spectra at all

    // careful: below may be end() at this point, yet below!=begin()
    if (below != begin()) --below; // we need to make one step left, so we are different from `above`
    // we are not at end() (or begin()==end())
    while (below != begin() && below->getMSLevel() != ms_level)
    {
      --below;
    }
    if (below->getMSLevel() != ms_level) return above; // below did not find anything valid; so it must be whatever `above` is (could be end())
    if (above == end()) return below;                  // queried beyond last spec, but we know there are spectra, so it must be whatever `below` is (which we know is valid)
    // we are between two spectra
    auto diff_left = RT - below->getRT();
    auto diff_right = above->getRT() - RT;
    return (diff_left < diff_right ? below : above);
  }

  MSExperiment::Iterator MSExperiment::getClosestSpectrumInRT(const double RT, UInt ms_level)
  {
    return begin() + std::distance(cbegin(), const_cast<const MSExperiment*>(this)->getClosestSpectrumInRT(RT, ms_level));
  }

  /// sets the chromatogram list
  void MSExperiment::setChromatograms(const std::vector<MSChromatogram > & chromatograms)
  {
    chromatograms_ = chromatograms;
  }

  /// sets the chromatogram list
  void MSExperiment::setChromatograms(std::vector<MSChromatogram> && chromatograms)
  {
    chromatograms_ = std::move(chromatograms);
  }

  /// adds a chromatogram to the list
  void MSExperiment::addChromatogram(const MSChromatogram & chromatogram)
  {
    chromatograms_.push_back(chromatogram);
  }

  void MSExperiment::addChromatogram(MSChromatogram&& chrom)
  {
    chromatograms_.push_back(std::move(chrom));
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
  const MSChromatogram MSExperiment::calculateTIC(float rt_bin_size, UInt ms_level) const
  {
    // The TIC is (re)calculated from the MS spectra with set ms_level (default 1).
    // Even if MSExperiment does not contain a TIC chromatogram explicitly, it can be reported.
    MSChromatogram TIC;
    for (const auto& spec: spectra_)
    {
      if ((spec.getMSLevel() == ms_level) || (ms_level == 0))
      {
        // fill chromatogram
        ChromatogramPeakType peak;
        peak.setRT(spec.getRT());
        peak.setIntensity(spec.calculateTIC());
        TIC.push_back(peak);
      }
    }
    if (rt_bin_size > 0)
    {
      LinearResamplerAlign lra;
      Param param = lra.getParameters();
      param.setValue("spacing", rt_bin_size);
      lra.setParameters(param);
      lra.raster(TIC);
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
    }
  }

  // static
  bool MSExperiment::containsScanOfLevel(size_t ms_level) const
  {
    //test if no scans with MS-level 1 exist
    for (const auto& spec : getSpectra())
    {
      if (spec.getMSLevel() == ms_level)
      {
        return true;
      }
    }
    return false;
  }

  bool MSExperiment::hasZeroIntensities(size_t ms_level) const
  {
    for (const auto& spec : getSpectra())
    {
      if (spec.getMSLevel() != ms_level)
      {
        continue;
      }
      for (const auto& p : spec)
      {
        if (p.getIntensity() == 0.0)
        {
          return true;
        }
      }
    }
    return false;
  }

  bool MSExperiment::hasPeptideIdentifications() const
  {
    for (const auto& spec : getSpectra())
    {
      if (!spec.getPeptideIdentifications().empty())
      {
        return true;
      }
    }
    return false;
  }

  bool MSExperiment::isIMFrame() const
  {
    if (spectra_.empty()) return false;
    auto rt_start = spectra_[0].getRT();
    auto last_drift = std::numeric_limits<double>::lowest();
    for (const auto& s : spectra_) {
      if (s.getRT() != rt_start) return false; // RT changes...
      if (s.getDriftTime() == last_drift) return false; // IM did not change...
      last_drift = s.getDriftTime();
    }
    return true; // RT stable, IM changing
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
  @brief Append a spectrum including float data arrays to current MSExperiment

  @param rt RT of new spectrum
  @param metadata_names Names of float data arrays attached to this spectrum
  @return Pointer to newly created spectrum
  */
  MSExperiment::SpectrumType* MSExperiment::createSpec_(PeakType::CoordinateType rt, const StringList& metadata_names)
  {
    SpectrumType* spectrum = createSpec_(rt);
    // create metadata arrays
    spectrum->getFloatDataArrays().reserve(metadata_names.size());
    for (StringList::const_iterator itm = metadata_names.begin(); itm != metadata_names.end(); ++itm)
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
    for (const MSSpectrum& spec : exp.getSpectra())
    {
      os << spec;
    }

    //chromatograms
    for (const MSChromatogram& chrom : exp.getChromatograms())
    {
      os << chrom;
    }

    os << "-- MSEXPERIMENT END --" << std::endl;

    return os;
  }
} //namespace OpenMS

