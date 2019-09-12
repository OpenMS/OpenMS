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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/AreaIterator.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <vector>


namespace OpenMS
{
  class Peak1D;
  class ChromatogramPeak;

  /**
    @brief In-Memory representation of a mass spectrometry experiment.

    Contains the data and metadata of an experiment performed with an MS (or
    HPLC and MS). This representation of an MS experiment is organized as list
    of spectra and chromatograms and provides an in-memory representation of
    popular mass-spectrometric file formats such as mzXML or mzML. The
    meta-data associated with an experiment is contained in
    ExperimentalSettings (by inheritance) while the raw data (as well as
    spectra and chromatogram level meta data) is stored in objects of type
    MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
    and getChromatogram functions.

    Be careful when changing the order of contained MSSpectrum instances, if
    tandem-MS data is stored in this class. The only way to find a precursor
    spectrum of MSSpectrum x is to search for the first spectrum before x that
    has a lower MS-level!

    @note For range operations, see \ref RangeUtils "RangeUtils module"!
    @note Some of the meta data is associated with the spectra directly (e.g. DataProcessing) and therefore the spectra need to be present to retain this information.
    @note For an on-disc representation of an MS experiment, see OnDiskExperiment.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MSExperiment :
    public RangeManager<2>,
    public ExperimentalSettings
  {

public:
    typedef Peak1D PeakT;
    typedef ChromatogramPeak ChromatogramPeakT;

    /// @name Base type definitions
    //@{
    /// Peak type
    typedef PeakT PeakType;
    /// Chromatogram peak type
    typedef ChromatogramPeakT ChromatogramPeakType;
    /// Area type
    typedef DRange<2> AreaType;
    /// Coordinate type of peak positions
    typedef PeakType::CoordinateType CoordinateType;
    /// Intensity type of peaks
    typedef PeakType::IntensityType IntensityType;
    /// RangeManager type
    typedef RangeManager<2> RangeManagerType;
    /// Spectrum Type
    typedef MSSpectrum SpectrumType;
    /// Chromatogram type
    typedef MSChromatogram ChromatogramType;
    /// STL base class type
    typedef std::vector<SpectrumType> Base;
    //@}

    /// @name Iterator type definitions
    //@{
    /// Mutable iterator
    typedef std::vector<SpectrumType>::iterator Iterator;
    /// Non-mutable iterator
    typedef std::vector<SpectrumType>::const_iterator ConstIterator;
    /// Mutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<PeakT, PeakT&, PeakT*, Iterator, SpectrumType::Iterator> AreaIterator;
    /// Immutable area iterator type (for traversal of a rectangular subset of the peaks)
    typedef Internal::AreaIterator<const PeakT, const PeakT&, const PeakT*, ConstIterator, SpectrumType::ConstIterator> ConstAreaIterator;
    //@}

    /// @name Delegations of calls to the vector of MSSpectra
    // Attention: these refer to the spectra vector only!
    //@{
    typedef Base::value_type value_type;
    typedef Base::iterator iterator;
    typedef Base::const_iterator const_iterator;

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
    void reserveSpaceSpectra(Size s);
    void reserveSpaceChromatograms(Size s);

    /// Constructor
    MSExperiment();

    /// Copy constructor
    MSExperiment(const MSExperiment & source);

    /// Move constructor
    MSExperiment(MSExperiment&&) = default;

    /// Assignment operator
    MSExperiment & operator=(const MSExperiment & source);

    /// Move assignment operator
    MSExperiment& operator=(MSExperiment&&) & = default;

    /// Assignment operator
    MSExperiment & operator=(const ExperimentalSettings & source);

    /// Equality operator
    bool operator==(const MSExperiment & rhs) const;

    /// Equality operator
    bool operator!=(const MSExperiment & rhs) const;

    ///@name Conversion to/from 2D data
    //@{
    /**
      @brief Reads out a 2D Spectrum

      Container can be a PeakArray or an STL container of peaks which
      supports push_back(), end() and back()
    */
    template <class Container>
    void get2DData(Container& cont) const
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
      @param store_metadata_names [MetaInfoInterface input only] Names of metadata arrays which should be created;
                                  data is filled from the metainfointerface of each element of the input container.
                                  Currently, only float data is supported!

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode)
    */
    template <class Container>
    void set2DData(const Container& container, const StringList& store_metadata_names)
    {
      // clean up the container first
      clear(true);
      SpectrumType* spectrum = nullptr;
      typename PeakType::CoordinateType current_rt = -std::numeric_limits<typename PeakType::CoordinateType>::max();
      for (typename Container::const_iterator iter = container.begin(); iter != container.end(); ++iter)
      {
        // check if the retention time has changed
        if (current_rt != iter->getRT() || spectrum == nullptr)
        {
          // append new spectrum
          if (current_rt > iter->getRT())
          {
            throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input container is not sorted!");
          }
          current_rt =  iter->getRT();
          spectrum = createSpec_(current_rt, store_metadata_names);
        }

        // add either data point or mass traces (depending on template argument value)
        ContainerAdd_<typename Container::value_type, false>::addData_(spectrum, &(*iter), store_metadata_names);
      }
    }

     /**
      @brief Assignment of a data container with RT and MZ to an MSExperiment

      Fill MSExperiment with data.
      Note that all data present (including meta-data) will be deleted prior to adding new data!

      @param container An iterable type whose elements support getRT(), getMZ() and getIntensity()
      @param add_mass_traces If true, each container element is searched for the metavalue
                             "num_of_masstraces".
                             If found, "masstrace_intensity" (X>=0) meta values are added as data points (with 13C spacing).
                             This is useful for, e.g., FF-Metabo output.
                             Note that the actual feature will NOT be added if mass traces are found (since MT0 is usually identical)

      @exception Exception::Precondition is thrown if the container is not sorted according to
      retention time (in debug AND release mode) OR a "masstrace_intensity" value is expected but not found
         
    */
    template <bool add_mass_traces, class Container>
    void set2DData(const Container& container)
    {
      // clean up the container first
      clear(true);
      SpectrumType* spectrum = nullptr;
      typename PeakType::CoordinateType current_rt = -std::numeric_limits<typename PeakType::CoordinateType>::max();
      for (typename Container::const_iterator iter = container.begin(); iter != container.end(); ++iter)
      {
        // check if the retention time has changed
        if (current_rt != iter->getRT() || spectrum == nullptr)
        {
          // append new spectrum
          if (current_rt > iter->getRT())
          {
            throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input container is not sorted!");
          }
          current_rt =  iter->getRT();
          spectrum = createSpec_(current_rt);
        }

        // add either data point or mass traces (depending on template argument value)
        ContainerAdd_<typename Container::value_type, add_mass_traces>::addData_(spectrum, &(*iter));
      }
    }

    //@}


    ///@name Iterating ranges and areas
    //@{
    /// Returns an area iterator for @p area
    AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz);

    /// Returns an invalid area iterator marking the end of an area
    AreaIterator areaEnd();

    /// Returns a non-mutable area iterator for @p area
    ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const;

    /// Returns an non-mutable invalid area iterator marking the end of an area
    ConstAreaIterator areaEndConst() const;

    /**
      @brief Fast search for spectrum range begin

      Returns the first scan which has equal or higher (>=) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTBegin(CoordinateType rt) const;

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      Returns the first scan which has higher (>) RT than @p rt.

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTEnd(CoordinateType rt) const;

    /**
      @brief Fast search for spectrum range begin

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTBegin(CoordinateType rt);

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      @note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTEnd(CoordinateType rt);

    //@}

    /**
      @name Range methods

      @note The range values (min, max, etc.) are not updated automatically. Call updateRanges() to update the values!
    */
    ///@{
    // Docu in base class
    void updateRanges() override;

    /**
      @brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level

      @param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
    */
    void updateRanges(Int ms_level);

    /// returns the minimal m/z value
    CoordinateType getMinMZ() const;

    /// returns the maximal m/z value
    CoordinateType getMaxMZ() const;

    /// returns the minimal retention time value
    CoordinateType getMinRT() const;

    /// returns the maximal retention time value
    CoordinateType getMaxRT() const;

    /**
      @brief Returns RT and m/z range the data lies in.

      RT is dimension 0, m/z is dimension 1
    */
    const AreaType& getDataRange() const;

    /// returns the total number of peaks
    UInt64 getSize() const;

    /// returns an array of MS levels
    const std::vector<UInt>& getMSLevels() const;

    ///@}

    ///@name Sorting spectra and peaks
    ///@{
    /**
      @brief Sorts the data points by retention time

      @param sort_mz if @em true, spectra are sorted by m/z position as well
    */
    void sortSpectra(bool sort_mz = true);

    /**
      @brief Sorts the data points of the chromatograms by m/z

      @param sort_rt if @em true, chromatograms are sorted by rt position as well
    */
    void sortChromatograms(bool sort_rt = true);

    /**
      @brief Checks if all spectra are sorted with respect to ascending RT

      @param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
    */
    bool isSorted(bool check_mz = true) const;

    //@}

    /// Resets all internal values
    void reset();

    /**
      @brief Clears the meta data arrays of all contained spectra (float, integer and string arrays)

      @return @em true if meta data arrays were present and removed. @em false otherwise.
    */
    bool clearMetaDataArrays();

    /// returns the meta information of this experiment (const access)
    const ExperimentalSettings& getExperimentalSettings() const;

    /// returns the meta information of this experiment (mutable access)
    ExperimentalSettings& getExperimentalSettings();

    /// get the file path to the first MS run
    void getPrimaryMSRunPath(StringList& toFill) const;

    /**
      @brief Returns the precursor spectrum of the scan pointed to by @p iterator

      If there is no precursor scan the past-the-end iterator is returned.
    */
    ConstIterator getPrecursorSpectrum(ConstIterator iterator) const;

    /// Swaps the content of this map with the content of @p from
    void swap(MSExperiment& from);

    /// sets the spectrum list
    void setSpectra(const std::vector<MSSpectrum>& spectra);

    /// adds a spectrum to the list
    void addSpectrum(const MSSpectrum& spectrum);

    void addSpectrum(MSSpectrum&& spectrum)
    {
      spectra_.push_back(std::forward<MSSpectrum>(spectrum));
    }

    /// returns the spectrum list
    const std::vector<MSSpectrum>& getSpectra() const;

    /// returns the spectrum list (mutable)
    std::vector<MSSpectrum>& getSpectra();

    /// sets the chromatogram list
    void setChromatograms(const std::vector<MSChromatogram>& chromatograms);

    /// adds a chromatogram to the list
    void addChromatogram(const MSChromatogram& chromatogram);

    void addChromatogram(MSChromatogram&& chrom)
    {
      chromatograms_.push_back(std::forward<MSChromatogram>(chrom));
    }

    /// returns the chromatogram list
    const std::vector<MSChromatogram>& getChromatograms() const;

    /// returns the chromatogram list (mutable)
    std::vector<MSChromatogram>& getChromatograms();

    /// @name Easy Access interface
    //@{
    /// returns a single chromatogram
    MSChromatogram& getChromatogram(Size id);

    /// returns a single spectrum
    MSSpectrum& getSpectrum(Size id);

    /// get the total number of spectra available
    Size getNrSpectra() const;

    /// get the total number of chromatograms available
    Size getNrChromatograms() const;
    //@}

    /**
    @brief Compute Total Ion Count per MS1 spectrum and applies the resampling algorithm, if a bin size in RT seconds greater than 0 is given.

    By default, each MS1 spectrum's intensity just gets summed up. Regular RT bins can be obtained by specifying @p rt_bin_size.

    @param bin_size RT bin size in seconds
    @return TIC Chromatogram
    **/
    const MSChromatogram getTIC(float rt_bin_size = 0) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

protected:

    /// MS levels of the data
    std::vector<UInt> ms_levels_;
    /// Number of all data points
    UInt64 total_size_;

    /// chromatograms
    std::vector<MSChromatogram > chromatograms_;

    /// spectra
    std::vector<SpectrumType> spectra_;

private:

    /// Helper class to add either general data points in set2DData or use mass traces from meta values
    template<typename ContainerValueType, bool addMassTraces>
    struct ContainerAdd_
    {
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item);
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item, const StringList& store_metadata_names);
    };

    template<typename ContainerValueType>
    struct ContainerAdd_<ContainerValueType, false>
    {
      /// general method for adding data points
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item)
      {
        // create temporary peak and insert it into spectrum
        spectrum->insert(spectrum->end(), PeakType());
        spectrum->back().setIntensity(item->getIntensity());
        spectrum->back().setPosition(item->getMZ());
      }
      /// general method for adding data points, including metadata arrays (populated from metainfointerface)
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item, const StringList& store_metadata_names)
      {
        addData_(spectrum, item);
        for (StringList::const_iterator itm = store_metadata_names.begin(); itm != store_metadata_names.end(); ++itm)
        {
          float val = std::numeric_limits<float>::quiet_NaN();
          if (item->metaValueExists(*itm)) val = item->getMetaValue(*itm);
          spectrum->getFloatDataArrays()[itm - store_metadata_names.begin()].push_back(val);
        }
      }
    };

    template<typename ContainerValueType>
    struct ContainerAdd_<ContainerValueType, true>
    {
      /// specialization for adding feature mass traces (does not support metadata_names currently)
      static void addData_(SpectrumType* spectrum, const ContainerValueType* item)
      {
        if (item->metaValueExists("num_of_masstraces"))
        {
          Size mts = item->getMetaValue("num_of_masstraces");
          int charge = (item->getCharge()==0 ? 1 : item->getCharge()); // set to 1 if charge is 0, otherwise div/0 below
          for (Size i = 0; i < mts; ++i)
          {
            String meta_name = String("masstrace_intensity_") + i;
            if (!item->metaValueExists(meta_name))
            {
              throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Meta value '") + meta_name + "' expected but not found in container.");
            }
            ContainerValueType p;
            p.setIntensity(item->getMetaValue(meta_name));
            p.setPosition(item->getMZ() + Constants::C13C12_MASSDIFF_U / charge * i);
            ContainerAdd_<ContainerValueType, false>::addData_(spectrum, &p);
          }
        }
        else ContainerAdd_<ContainerValueType, false>::addData_(spectrum, item);
      }
    };

    /*
      @brief Append a spectrum to current MSExperiment

      @param rt RT of new spectrum
      @return Pointer to newly created spectrum
    */
    SpectrumType* createSpec_(PeakType::CoordinateType rt);

    /*
      @brief Append a spectrum including floatdata arrays to current MSExperiment

      @param rt RT of new spectrum
      @param metadata_names Names of floatdata arrays attached to this spectrum
      @return Pointer to newly created spectrum
    */
    SpectrumType* createSpec_(PeakType::CoordinateType rt, const StringList& metadata_names);

  };

  /// Print the contents to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const MSExperiment& exp);

} // namespace OpenMS

#include <OpenMS/KERNEL/StandardTypes.h>


