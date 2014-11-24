// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSEXPERIMENT_H
#define OPENMS_KERNEL_MSEXPERIMENT_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/AreaIterator.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <vector>
#include <algorithm>
#include <limits>

namespace OpenMS
{
  class Peak1D;

  /**
    @brief In-Memory representation of a mass spectrometry experiment.

    Contains the data and metadata of an experiment performed with an MS (or HPLC and MS). This representation of an MS experiment is organized as list of spectra and chromatograms and provides an in-memory representation of popular mass-spectrometric file formats such as mzXML or mzML. The meta-data associated with an experiment is contained in ExperimentalSettings (by inheritance) while the raw data (as well as spectra and chromatogram level meta data) is stored in objects of type MSSpectrum and MSChromatogram which are accessible through the getSpectrum and getChromatogam functions.

    Be careful when changing the order of contained MSSpectrum instances, if tandem-MS data is
    stored in this class. The only way to find a precursor spectrum of MSSpectrum x is to
    search for the first spectrum before x that has a lower MS-level!

    @note For range operations, see \ref RangeUtils "RangeUtils module"!
    @note Some of the meta data is associated with the spectra directly (e.g. DataProcessing) and therefore the spectra need to be present to retain this information.
    @note For an on-disc representation of an MS experiment, see OnDiskExperiment.

    @ingroup Kernel
  */
  template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
  class MSExperiment :
    public RangeManager<2>,
    public ExperimentalSettings
  {

public:

    /// @name Base type definitions
    //@{
    /// Peak type
    typedef PeakT PeakType;
    /// Chromatogram peak type
    typedef ChromatogramPeakT ChromatogramPeakType;
    /// Area type
    typedef DRange<2> AreaType;
    /// Coordinate type of peak positions
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Intensity type of peaks
    typedef typename PeakType::IntensityType IntensityType;
    /// RangeManager type
    typedef RangeManager<2> RangeManagerType;
    /// Spectrum Type
    typedef MSSpectrum<PeakType> SpectrumType;
    /// Chromatogram type
    typedef MSChromatogram<ChromatogramPeakType> ChromatogramType;
    /// STL base class type
    typedef std::vector<SpectrumType> Base;
    //@}

    /// @name Iterator type definitions
    //@{
    /// Mutable iterator
    typedef typename std::vector<SpectrumType>::iterator Iterator;
    /// Non-mutable iterator
    typedef typename std::vector<SpectrumType>::const_iterator ConstIterator;
    /// Mutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<PeakT, PeakT &, PeakT *, Iterator, typename SpectrumType::Iterator> AreaIterator;
    /// Immutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<const PeakT, const PeakT &, const PeakT *, ConstIterator, typename SpectrumType::ConstIterator> ConstAreaIterator;
    //@}

    /// @name Delegations of calls to the vector of MSSpectra
    // Attention: these refer to the spectra vector only!
    //@{
    typedef typename Base::value_type value_type; 
    typedef typename Base::iterator iterator; 
    typedef typename Base::const_iterator const_iterator; 

    inline Size size() const
    {
      return spectra_.size(); 
    }

    inline void resize(Size s)
    {
      spectra_.resize(s); 
    }

    inline bool empty() const
    {
      return spectra_.empty(); 
    }

    inline void reserve(Size s)
    {
      spectra_.reserve(s); 
    }

    inline SpectrumType& operator[] (Size n)
    {
      return spectra_[n];
    }

    inline const SpectrumType& operator[] (Size n) const
    {
      return spectra_[n];
    }

    inline Iterator begin() 
    {
      return spectra_.begin();
    }

    inline ConstIterator begin() const
    {
      return spectra_.begin();
    }

    inline Iterator end() 
    {
      return spectra_.end();
    }

    inline ConstIterator end() const
    {
      return spectra_.end();
    }
    //@}

    // Aliases / chromatograms
    inline void reserveSpaceSpectra(Size s)
    {
      spectra_.reserve(s); 
    }
    inline void reserveSpaceChromatograms(Size s)
    {
      chromatograms_.reserve(s); 
    }

    /// Constructor
    MSExperiment() :
      RangeManagerType(),
      ExperimentalSettings(),
      ms_levels_(),
      total_size_(0)
    {}

    /// Copy constructor
    MSExperiment(const MSExperiment & source) :
      RangeManagerType(source),
      ExperimentalSettings(source),
      ms_levels_(source.ms_levels_),
      total_size_(source.total_size_),
      chromatograms_(source.chromatograms_),
      spectra_(source.spectra_)
    {}

    /// Assignment operator
    MSExperiment & operator=(const MSExperiment & source)
    {
      if (&source == this) return *this;

      RangeManagerType::operator=(source);
      ExperimentalSettings::operator=(source);

      ms_levels_     = source.ms_levels_;
      total_size_    = source.total_size_;
      chromatograms_ = source.chromatograms_;
      spectra_ = source.spectra_;

      //no need to copy the alloc?!
      //alloc_

      return *this;
    }

    /// Assignment operator
    MSExperiment & operator=(const ExperimentalSettings & source)
    {
      ExperimentalSettings::operator=(source);
      return *this;
    }

    /// Equality operator
    bool operator==(const MSExperiment & rhs) const
    {
      return ExperimentalSettings::operator==(rhs) &&
          chromatograms_ == rhs.chromatograms_ && 
          spectra_ == rhs.spectra_;
    }

    /// Equality operator
    bool operator!=(const MSExperiment & rhs) const
    {
      return !(operator==(rhs));
    }

    ///@name Conversion to/from 2D data
    //@{
    /**
      @brief Reads out a 2D Spectrum

      Container can be a PeakArray or an STL container of peaks which
      supports push_back(), end() and back()
    */
    template <class Container>
    void get2DData(Container & cont) const
    {
      for (typename Base::const_iterator spec = spectra_.begin(); spec != spectra_.end(); ++spec)
      {
        if (spec->getMSLevel() != 1)
        {
          continue;
        }
				typename Container::value_type s; // explicit object here, since instantiation within push_back() fails on VS<12
        for (typename SpectrumType::const_iterator it = spec->begin(); it != spec->end(); ++it)
        {
          cont.push_back(s);
          cont.back().setRT(spec->getRT());
          cont.back().setMZ(it->getMZ());
          cont.back().setIntensity(it->getIntensity());
        }
      }
    }

    /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @param container An iterable type whose elements support getRT(), getMZ() and getIntensity()

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode)
    */
    template <class Container>
    void set2DData(const Container& container)
    {
      set2DData<false, Container>(container);
    }

     /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @param container An iterable type whose elements support getRT(), getMZ() and getIntensity()
      @param add_mass_traces If true, each container element is searched for the metavalue
                             "num_of_masstraces".
                             If found, "masstrace_intensity_<X>" (X>=0) meta values are added as data points (with 13C spacing).
                             This is useful for, e.g., FF-Metabo output.
                             Note that the actual feature will NOT be added if mass traces are found (since MT0 is usually identical)

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode) OR a "masstrace_intensity_<X>" value is expected but not found
         
    */
    template <bool add_mass_traces, class Container>
    void set2DData(const Container& container)
    {
      // clean up the container first
      clear(true);

      SpectrumType* spectrum = 0;

      typename PeakType::CoordinateType current_rt = -std::numeric_limits<typename PeakType::CoordinateType>::max();

      for (typename Container::const_iterator iter = container.begin(); iter != container.end(); ++iter)
      {
        // check if the retention time has changed
        if (current_rt != iter->getRT() || spectrum == 0)
        {
          // append new spectrum
          if (current_rt > iter->getRT())
          {
            throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input container is not sorted!");
          }
          current_rt =  iter->getRT();
          spectra_.insert(spectra_.end(), SpectrumType());
          spectrum = &(spectra_.back());
          spectrum->setRT(current_rt);
          spectrum->setMSLevel(1);
        }

        // add either datapoint or mass traces (depending on template argument value)
        ContainerAdd_<typename Container::const_iterator, add_mass_traces>::addData_(spectrum, iter);
      }
    }

    //@}


    ///@name Iterating ranges and areas
    //@{
    /// Returns an area iterator for @p area
    AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz)
    {
      OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
      OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
      OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using AreaIterator will give invalid results!")
      //std::cout << "areaBegin: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
      return AreaIterator(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
    }

    /// Returns an invalid area iterator marking the end of an area
    AreaIterator areaEnd()
    {
      return AreaIterator();
    }

    /// Returns a non-mutable area iterator for @p area
    ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const
    {
      OPENMS_PRECONDITION(min_rt <= max_rt, "Swapped RT range boundaries!")
      OPENMS_PRECONDITION(min_mz <= max_mz, "Swapped MZ range boundaries!")
      OPENMS_PRECONDITION(this->isSorted(true), "Experiment is not sorted by RT and m/z! Using ConstAreaIterator will give invalid results!")
      //std::cout << "areaBeginConst: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
      return ConstAreaIterator(spectra_.begin(), RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
    }

    /// Returns an non-mutable invalid area iterator marking the end of an area
    ConstAreaIterator areaEndConst() const
    {
      return ConstAreaIterator();
    }

    /**
      @brief Fast search for spectrum range begin

      Returns the first scan which has equal or higher (>=) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTBegin(CoordinateType rt) const
    {
      SpectrumType s;
      s.setRT(rt);
      return lower_bound(spectra_.begin(), spectra_.end(), s, typename SpectrumType::RTLess());
    }

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      Returns the first scan which has higher (>) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTEnd(CoordinateType rt) const
    {
      SpectrumType s;
      s.setRT(rt);
      return upper_bound(spectra_.begin(), spectra_.end(), s, typename SpectrumType::RTLess());
    }

    /**
      @brief Fast search for spectrum range begin

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTBegin(CoordinateType rt)
    {
      SpectrumType s;
      s.setRT(rt);
      return lower_bound(spectra_.begin(), spectra_.end(), s, typename SpectrumType::RTLess());
    }

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTEnd(CoordinateType rt)
    {
      SpectrumType s;
      s.setRT(rt);
      return upper_bound(spectra_.begin(), spectra_.end(), s, typename SpectrumType::RTLess());
    }

    //@}

    /**
      @name Range methods

      @note The range values (min, max, etc.) are not updated automatically. Call updateRanges() to update the values!
    */
    ///@{
    // Docu in base class
    virtual void updateRanges()
    {
      updateRanges(-1);
    }

    /**
      @brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level

      @param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
    */
    void updateRanges(Int ms_level)
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
      for (typename Base::iterator it = spectra_.begin(); it != spectra_.end(); ++it)
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

      for (typename std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
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
    CoordinateType getMinMZ() const
    {
      return RangeManagerType::pos_range_.minPosition()[1];
    }

    /// returns the maximal m/z value
    CoordinateType getMaxMZ() const
    {
      return RangeManagerType::pos_range_.maxPosition()[1];
    }

    /// returns the minimal retention time value
    CoordinateType getMinRT() const
    {
      return RangeManagerType::pos_range_.minPosition()[0];
    }

    /// returns the maximal retention time value
    CoordinateType getMaxRT() const
    {
      return RangeManagerType::pos_range_.maxPosition()[0];
    }

    /**
      @brief Returns RT and m/z range the data lies in.

      RT is dimension 0, m/z is dimension 1
    */
    const AreaType & getDataRange() const
    {
      return RangeManagerType::pos_range_;
    }

    /// returns the total number of peaks
    UInt64 getSize() const
    {
      return total_size_;
    }

    /// returns an array of MS levels
    const std::vector<UInt> & getMSLevels() const
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
    void sortSpectra(bool sort_mz = true)
    {
      std::sort(spectra_.begin(), spectra_.end(), typename SpectrumType::RTLess());

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
    void sortChromatograms(bool sort_rt = true)
    {
      // sort the chromatograms according to their product m/z
      std::sort(chromatograms_.begin(), chromatograms_.end(), typename ChromatogramType::MZLess());

      if (sort_rt)
      {
        for (typename std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
        {
          it->sortByPosition();
        }
      }
    }

    /**
      @brief Checks if all spectra are sorted with respect to ascending RT

      @param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
    */
    bool isSorted(bool check_mz = true) const
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
    void reset()
    {
      spectra_.clear();           //remove data
      RangeManagerType::clearRanges();           //reset range manager
      ExperimentalSettings::operator=(ExperimentalSettings());           //reset meta info
    }

    /**
      @brief Clears the meta data arrays of all contained spectra (float, integer and string arrays)

      @return @em true if meta data arrays were present and removed. @em false otherwise.
    */
    bool clearMetaDataArrays()
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
    const ExperimentalSettings & getExperimentalSettings() const
    {
      return *this;
    }

    /// returns the meta information of this experiment (mutable access)
    ExperimentalSettings & getExperimentalSettings()
    {
      return *this;
    }

    /**
      @brief Returns the precursor spectrum of the scan pointed to by @p iterator

      If there is no precursor scan the past-the-end iterator is returned.
    */
    ConstIterator getPrecursorSpectrum(ConstIterator iterator) const
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
      }
      while (iterator != spectra_.begin());

      return spectra_.end();
    }

    /// Swaps the content of this map with the content of @p from
    void swap(MSExperiment & from)
    {
      MSExperiment tmp;

      //swap range information
      tmp.RangeManagerType::operator=(* this);
      this->RangeManagerType::operator=(from);
      from.RangeManagerType::operator=(tmp);

      //swap experimental settings
      tmp.ExperimentalSettings::operator=(* this);
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

    /// sets the spectra list
    void setSpectra(const std::vector<MSSpectrum<PeakT> > & spectra)
    {
      spectra_ = spectra;
    }

    /// adds a spectra to the list
    void addSpectrum(const MSSpectrum<PeakT> & spectrum)
    {
      spectra_.push_back(spectrum);
    }

    /// returns the spectra list
    const std::vector<MSSpectrum<PeakT> > & getSpectra() const
    {
      return spectra_;
    }

    /// returns the spectra list
    std::vector<MSSpectrum<PeakT> > & getSpectra() 
    {
      return spectra_;
    }

    /// sets the chromatogram list
    void setChromatograms(const std::vector<MSChromatogram<ChromatogramPeakType> > & chromatograms)
    {
      chromatograms_ = chromatograms;
    }

    /// adds a chromatogram to the list
    void addChromatogram(const MSChromatogram<ChromatogramPeakType> & chromatogram)
    {
      chromatograms_.push_back(chromatogram);
    }

    /// returns the chromatogram list
    const std::vector<MSChromatogram<ChromatogramPeakType> > & getChromatograms() const
    {
      return chromatograms_;
    }

    /// @name Easy Access interface
    //@{
    /// returns a single chromatogram 
    MSChromatogram<ChromatogramPeakType> & getChromatogram(Size id) 
    {
      return chromatograms_[id];
    }

    /// returns a single spectrum 
    MSSpectrum<PeakT> & getSpectrum(Size id) 
    {
      return spectra_[id];
    }

    /// get the total number of spectra available
    inline Size getNrSpectra() const
    {
      return spectra_.size();
    }

    /// get the total number of chromatograms available
    inline Size getNrChromatograms() const
    {
      return chromatograms_.size();
    }
    //@}

    /// returns the total ion chromatogram (TIC)
    const MSChromatogram<ChromatogramPeakType> getTIC() const
    {
      // The TIC is (re)calculated from the MS1 spectra. Even if MSExperiment does not contain a TIC chromatogram explicitly, it can be reported.
      MSChromatogram<ChromatogramPeakType> TIC;
      for (typename Base::const_iterator spec_it = spectra_.begin(); spec_it != spectra_.end(); ++spec_it)
      {
        if (spec_it->getMSLevel() == 1)
        {
          double totalIntensity = 0;
          // sum intensities of a spectrum
          for (typename SpectrumType::const_iterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
          {
            totalIntensity += peak_it->getIntensity();
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
    void clear(bool clear_meta_data)
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

protected:


    /// MS levels of the data
    std::vector<UInt> ms_levels_;
    /// Number of all data points
    UInt64 total_size_;

    /// chromatograms
    std::vector<MSChromatogram<ChromatogramPeakType> > chromatograms_;

    /// spectra
    std::vector<SpectrumType> spectra_;

private:
   
    /// Helper class to add either general data points in set2DData or use mass traces from meta values
    template<typename ContainerIterator, bool addMassTraces>
    struct ContainerAdd_
    {
      static void addData_(SpectrumType* spectrum, ContainerIterator& iter);      
    };

    template<typename ContainerIterator>
    struct ContainerAdd_<ContainerIterator, false>
    {
      /// general method for adding data points (no mass traces desired or found)
      static void addData_(SpectrumType* spectrum, ContainerIterator& iter)
      {
        // create temporary peak and insert it into spectrum
        spectrum->insert(spectrum->end(), PeakType());
        spectrum->back().setIntensity(iter->getIntensity());
        spectrum->back().setPosition(iter->getMZ());
      }
    };

    template<typename ContainerIterator>
    struct ContainerAdd_<ContainerIterator, true>
    {
      /// specialization for adding feature mass traces
      static void addData_(SpectrumType* spectrum, ContainerIterator& iter)
      {
        if (iter->metaValueExists("num_of_masstraces"))
        {
          Size mts = iter->getMetaValue("num_of_masstraces");
          for (Size i=0; i<mts; ++i)
          {
            String meta_name = String("masstrace_intensity_") + i;
            if (!iter->metaValueExists(meta_name))
            {
              throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Meta value '") + meta_name + "' expected but not found in container.");
            }
            spectrum->insert(spectrum->end(), PeakType());
            spectrum->back().setIntensity(iter->getMetaValue(meta_name));
            spectrum->back().setPosition(iter->getMZ() + Constants::C13C12_MASSDIFF_U/iter->getCharge()*i);
          }
        }
        else ContainerAdd_<ContainerIterator, false>::addData_(spectrum, iter);
      }
    };

  };


  /// Print the contents to a stream.
  template <typename PeakT, typename ChromatogramPeakT>
  std::ostream & operator<<(std::ostream & os, const MSExperiment<PeakT, ChromatogramPeakT> & exp)
  {
    os << "-- MSEXPERIMENT BEGIN --" << std::endl;

    //experimental settings
    os << static_cast<const ExperimentalSettings &>(exp);

    //spectra
    for (typename MSExperiment<PeakT>::const_iterator it = exp.begin(); it != exp.end(); ++it)
    {
      os << *it;
    }

    //chromatograms
    for (typename std::vector<MSChromatogram<ChromatogramPeakT> >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
    {
      os << *it;
    }

    os << "-- MSEXPERIMENT END --" << std::endl;

    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSEXPERIMENT_H
