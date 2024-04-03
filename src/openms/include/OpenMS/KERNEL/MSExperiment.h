// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Exception.h>
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
    @brief In-Memory representation of a mass spectrometry run.

    This representation of an MS run is organized as list
    of spectra and chromatograms and provides an in-memory representation of
    popular mass-spectrometric file formats such as mzXML or mzML. The
    meta-data associated with an experiment is contained in
    ExperimentalSettings (by inheritance) while the raw data (as well as
    spectra and chromatogram level meta data) is stored in objects of type
    MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
    and getChromatogram functions.

    @note For range operations, see \ref RangeUtils "RangeUtils module"!
    @note Some of the meta data is associated with the spectra directly (e.g. DataProcessing) and therefore the spectra need to be present to retain this information.
    @note For an on-disc representation of an MS experiment, see OnDiskExperiment.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MSExperiment final : public RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility>,
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
    /// Coordinate type of peak positions
    typedef PeakType::CoordinateType CoordinateType;
    /// Intensity type of peaks
    typedef PeakType::IntensityType IntensityType;
    /// RangeManager type
    typedef RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility> RangeManagerType;
    /// RangeManager type
    typedef RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility> RangeManagerContainerType;
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

    /// D'tor
    ~MSExperiment() override;

    /// Equality operator
    bool operator==(const MSExperiment & rhs) const;

    /// Equality operator
    bool operator!=(const MSExperiment & rhs) const;
    
    /// The number of spectra
    inline Size size() const noexcept
    {
      return spectra_.size();
    }

    /// Resize to @p n spectra
    inline void resize(Size n)
    {
      spectra_.resize(n);
    }

    /// Are there any spectra (does not consider chromatograms)
    inline bool empty() const noexcept
    {
      return spectra_.empty();
    }
    
    /// Reserve space for @p n spectra
    inline void reserve(Size n)
    {
      spectra_.reserve(n);
    }

    /// Random access to @p n'th spectrum
    inline SpectrumType& operator[](Size n)
    {
      return spectra_[n];
    }

    /// Random access to @p n'th spectrum
    inline const SpectrumType& operator[](Size n) const
    {
      return spectra_[n];
    }

    inline Iterator begin() noexcept
    {
      return spectra_.begin();
    }

    inline ConstIterator begin() const noexcept
    {
      return spectra_.cbegin();
    }

    inline ConstIterator cbegin() const noexcept
    {
      return spectra_.cbegin();
    }

    inline Iterator end()
    {
      return spectra_.end();
    }

    inline ConstIterator end() const noexcept
    {
      return spectra_.cend();
    }
    
    inline ConstIterator cend() const noexcept
    {
      return spectra_.cend();
    }
    //@}

    // Aliases / chromatograms
    void reserveSpaceSpectra(Size s);
    void reserveSpaceChromatograms(Size s);

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

      @tparam Container An iterable type whose elements support getRT(), getMZ() and getIntensity()
      @tparam add_mass_traces If true, each container element is searched for the metavalue
                             "num_of_masstraces".
                             If found, "masstrace_intensity" (X>=0) meta values are added as data points (with 13C spacing).
                             This is useful for, e.g., FF-Metabo output.
                             Note that the actual feature will NOT be added if mass traces are found (since MT0 is usually identical)
      @param container The input data with RT,m/z and intensity

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
    AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt  ms_level = 1);

    /// Returns a area iterator for all peaks in @p range. If a dimension is empty(), it is ignored (i.e. does not restrict the area)
    AreaIterator areaBegin(const RangeManagerType& range, UInt ms_level = 1);

    /// Returns an invalid area iterator marking the end of an area
    AreaIterator areaEnd();

    /// Returns a non-mutable area iterator for @p area
    ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt ms_level = 1) const;

    /// Returns a non-mutable area iterator for all peaks in @p range. If a dimension is empty(), it is ignored (i.e. does not restrict the area)
    ConstAreaIterator areaBeginConst(const RangeManagerType& range, UInt ms_level = 1) const;

    /// Returns an non-mutable invalid area iterator marking the end of an area
    ConstAreaIterator areaEndConst() const;

    // for fast pyOpenMS access to MS1 peak data in format: [rt, [mz, intensity]]
    void get2DPeakData(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, 
      std::vector<float>& rt, 
      std::vector<std::vector<float>>& mz, 
      std::vector<std::vector<float>>& intensity) const
    {
      float t = -1.0;
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz); it != areaEndConst(); ++it)
      {
        if (it.getRT() != t) 
        {
          t = (float)it.getRT();
          rt.push_back(t);
        }
        mz.back().push_back((float)it->getMZ());
        intensity.back().push_back(it->getIntensity());
      }
    }

    // for fast pyOpenMS access to MS1 peak data in format: [rt, [mz, intensity, ion mobility]]
    void get2DPeakData(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, 
      std::vector<float>& rt, 
      std::vector<std::vector<float>>& mz,
      std::vector<std::vector<float>>& intensity, 
      std::vector<std::vector<float>>& ion_mobility) const
    {
      float t = -1.0;
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz); it != areaEndConst(); ++it)
      {
        if (it.getRT() != t)
        {
          t = (float)it.getRT();
          rt.push_back(t);
        }
        
        const MSSpectrum& spectrum = it.getSpectrum();
        bool has_IM = spectrum.containsIMData();
        float peak_IM{-1.0f};
        if (has_IM)
        {
          const auto& im_data = spectrum.getIMData();
          const Size peak_index = it.getPeakIndex().peak;
          if (spectrum.getFloatDataArrays()[im_data.first].size() == spectrum.size())
          {
            peak_IM = spectrum.getFloatDataArrays()[im_data.first][peak_index];
          }          
        }
        ion_mobility.back().push_back(peak_IM);
        mz.back().push_back((float)it->getMZ());
        intensity.back().push_back(it->getIntensity());
      }
    }

    // for fast pyOpenMS access to MS1 peak data in format: [rt, mz, intensity]
    void get2DPeakData(
      CoordinateType min_rt,
      CoordinateType max_rt,
      CoordinateType min_mz,
      CoordinateType max_mz,
      std::vector<float>& rt,
      std::vector<float>& mz,
      std::vector<float>& intensity) const
    {
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz); it != areaEndConst(); ++it)
      {
        rt.push_back((float)it.getRT());
        mz.push_back((float)it->getMZ());
        intensity.push_back(it->getIntensity());
      }
    }

    // for fast pyOpenMS access to MS1 peak data in format: [rt, mz, intensity, ion mobility]
    void get2DPeakDataIon(
      CoordinateType min_rt,
      CoordinateType max_rt,
      CoordinateType min_mz,
      CoordinateType max_mz,
      std::vector<float>& rt,
      std::vector<float>& mz,
      std::vector<float>& intensity,
      std::vector<float>& ion_mobility) const
    {
      for (auto it = areaBeginConst(min_rt, max_rt, min_mz, max_mz); it != areaEndConst(); ++it)
      {
        rt.push_back((float)it.getRT());
        mz.push_back((float)it->getMZ());
        intensity.push_back(it->getIntensity());

        const MSSpectrum& spectrum = it.getSpectrum();
        bool has_IM = spectrum.containsIMData();
        float peak_IM = -1.0;
        if (has_IM)
        {
          const auto& im_data = spectrum.getIMData();
          const Size& peak_index = it.getPeakIndex().peak;
          if (spectrum.getFloatDataArrays()[im_data.first].size() == spectrum.size())
          {
            peak_IM = spectrum.getFloatDataArrays()[im_data.first][peak_index];
          }          
        }        
        ion_mobility.push_back(peak_IM);
      }
    }

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


    /**
      @brief Fast search for spectrum range begin

      Returns the first scan which has equal or higher (>=) ion mobility than @p rt.

      @note Make sure the spectra are sorted with respect to ion mobility! Otherwise the result is undefined.
    */
    ConstIterator IMBegin(CoordinateType im) const;

    /**
      @brief Fast search for spectrum range end (returns the past-the-end iterator)

      Returns the first scan which has higher (>) ion mobility than @p im.

      @note Make sure the spectra are sorted with respect to ion mobility! Otherwise the result is undefined.
    */
    ConstIterator IMEnd(CoordinateType im) const;
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

    /// returns the total number of peaks
    UInt64 getSize() const;

    /// returns an array of MS levels
    const std::vector<UInt>& getMSLevels() const;

    ///@}

    /// If the file is loaded from an sqMass file, this run-ID allows to connect to the corresponding OSW identification file
    /// If the run-ID was not stored (older version) or this MSExperiment was not loaded from sqMass, then 0 is returned.
    UInt64 getSqlRunID() const;

    /// sets the run-ID which is used when storing an sqMass file
    void setSqlRunID(UInt64 id);

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

    /// Clear all internal data (spectra, ranges, metadata)
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

    /**
      @brief Returns the index of the precursor spectrum for spectrum at index @p zero_based_index

      If there is no precursor scan -1 is returned.
    */
    int getPrecursorSpectrum(int zero_based_index) const;

    /// Swaps the content of this map with the content of @p from
    void swap(MSExperiment& from);

    /// sets the spectrum list
    void setSpectra(const std::vector<MSSpectrum>& spectra);
    void setSpectra(std::vector<MSSpectrum>&& spectra);

    /// adds a spectrum to the list
    void addSpectrum(const MSSpectrum& spectrum);
    void addSpectrum(MSSpectrum&& spectrum);

    /// returns the spectrum list
    const std::vector<MSSpectrum>& getSpectra() const;

    /// returns the spectrum list (mutable)
    std::vector<MSSpectrum>& getSpectra();

    /// Returns the closest(=nearest) spectrum in retention time to the given RT
    ConstIterator getClosestSpectrumInRT(const double RT) const;
    Iterator getClosestSpectrumInRT(const double RT);

    /// Returns the closest(=nearest) spectrum in retention time to the given RT of a certain MS level
    ConstIterator getClosestSpectrumInRT(const double RT, UInt ms_level) const;
    Iterator getClosestSpectrumInRT(const double RT, UInt ms_level);

    /// sets the chromatogram list
    void setChromatograms(const std::vector<MSChromatogram>& chromatograms);
    void setChromatograms(std::vector<MSChromatogram>&& chromatograms);

    /// adds a chromatogram to the list
    void addChromatogram(const MSChromatogram& chromatogram);
    void addChromatogram(MSChromatogram&& chrom);

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
    @brief Computes the total ion chromatogram (TIC) for a given MS level (use ms_level = 0 for all levels). 

    By default, each MS spectrum's intensity just gets summed up. Regular RT bins can be obtained by specifying @p rt_bin_size.
    If a bin size in RT seconds greater than 0 is given resampling is used.

    @param rt_bin_size RT bin size in seconds (0 = no resampling)
    @param ms_level MS level of spectra for calculation (0 = all levels)
    @return TIC Chromatogram
    **/
    const MSChromatogram calculateTIC(float rt_bin_size = 0, UInt ms_level = 1) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    /// returns true if at least one of the spectra has the specified level
    bool containsScanOfLevel(size_t ms_level) const;

    /// returns true if any MS spectra of trthe specified level contain at least one peak with intensity of 0.0
    bool hasZeroIntensities(size_t ms_level) const;

    /// do any of the spectra have a PeptideID?
    bool hasPeptideIdentifications() const;

    /// Are all MSSpectra in this experiment part of an IM Frame? I.e. they all have the same RT, but different drift times
    bool isIMFrame() const;

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


