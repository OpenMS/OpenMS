// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <numeric>

namespace OpenMS
{
  enum class DriftTimeUnit;
  /**
    @brief The representation of a 1D spectrum.

    It contains peak data and metadata about specific instrument settings,
    acquisition settings, description of the meta values used in the peaks and precursor info
    (SpectrumSettings).

    Several MSSpectrum instances are contained in a peak map (MSExperiment), which is essentially
    a vector of spectra with additional information about the experiment.

    Precursor info from SpectrumSettings should only be used if this spectrum is a tandem-MS
    spectrum. The precursor spectrum is the first spectrum in MSExperiment, that has a lower
    MS-level than the current spectrum.

    @note For range operations, see \ref RangeUtils "RangeUtils module"!

    @ingroup Kernel
  */
  class OPENMS_DLLAPI MSSpectrum final :
    private std::vector<Peak1D>,
    public RangeManagerContainer<RangeMZ, RangeIntensity, RangeMobility>,
    public SpectrumSettings
  {
public:

    /// Comparator for the retention time.
    struct OPENMS_DLLAPI RTLess
    {
      bool operator()(const MSSpectrum& a, const MSSpectrum& b) const;
    };
    /// Comparator for the ion mobility.
    struct OPENMS_DLLAPI IMLess {
      bool operator()(const MSSpectrum& a, const MSSpectrum& b) const;
    };

    /// Used to remember what subsets in a spectrum are sorted already to allow faster sorting of the spectrum
    struct Chunk {
      Size start; ///< inclusive
      Size end; ///< not inclusive
      bool is_sorted; ///< are the Peaks in [start, end) sorted yet?
      Chunk(Size p_start, Size p_end, bool p_sorted) : start(p_start), end(p_end), is_sorted(p_sorted)
      {
      }
    };

    struct Chunks {
      public:
        Chunks(const MSSpectrum& s) : spec_(s) {}
        void add(bool is_sorted)
        {
          chunks_.emplace_back((chunks_.empty() ? 0 : chunks_.back().end), spec_.size(), is_sorted);
        }
        std::vector<Chunk>& getChunks()
        {
          return chunks_;
        }
      private:
        std::vector<Chunk> chunks_;
        const MSSpectrum& spec_;
    };

    ///@name Base type definitions
    //@{
    /// Peak type
    typedef OpenMS::Peak1D PeakType;
    /// Coordinate (m/z) type
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Spectrum base type
    typedef std::vector<PeakType> ContainerType;
    /// RangeManager
    typedef RangeManagerContainer<RangeMZ, RangeIntensity, RangeMobility> RangeManagerContainerType;
    typedef RangeManager<RangeMZ, RangeIntensity, RangeMobility> RangeManagerType;
    /// Float data array vector type
    typedef OpenMS::DataArrays::FloatDataArray FloatDataArray ;
    typedef std::vector<FloatDataArray> FloatDataArrays;
    /// String data array vector type
    typedef OpenMS::DataArrays::StringDataArray StringDataArray ;
    typedef std::vector<StringDataArray> StringDataArrays;
    /// Integer data array vector type
    typedef OpenMS::DataArrays::IntegerDataArray IntegerDataArray ;
    typedef std::vector<IntegerDataArray> IntegerDataArrays;
    //@}

    ///@name Peak container iterator type definitions
    //@{
    /// Mutable iterator
    typedef typename ContainerType::iterator Iterator;
    /// Non-mutable iterator
    typedef typename ContainerType::const_iterator ConstIterator;
    /// Mutable reverse iterator
    typedef typename ContainerType::reverse_iterator ReverseIterator;
    /// Non-mutable reverse iterator
    typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
    //@}

    ///@name Export methods from std::vector<Peak1D>
    //@{
    using ContainerType::operator[];
    using ContainerType::begin;
    using ContainerType::rbegin;
    using ContainerType::end;
    using ContainerType::rend;
    using ContainerType::cbegin;
    using ContainerType::cend;
    using ContainerType::resize;
    using ContainerType::size;
    using ContainerType::push_back;
    using ContainerType::emplace_back;
    using ContainerType::pop_back;
    using ContainerType::empty;
    using ContainerType::front;
    using ContainerType::back;
    using ContainerType::reserve;
    using ContainerType::insert;
    using ContainerType::erase;
    using ContainerType::swap;
    using ContainerType::data;
    using ContainerType::shrink_to_fit;

    using typename ContainerType::iterator;
    using typename ContainerType::const_iterator;
    using typename ContainerType::size_type;
    using typename ContainerType::value_type;
    using typename ContainerType::reference;
    using typename ContainerType::const_reference;
    using typename ContainerType::pointer;
    using typename ContainerType::difference_type;

    //@}


    /// Constructor
    MSSpectrum();

    /// Constructor from a list of Peak1D, e.g. MSSpectrum spec{ {mz1, int1}, {mz2, int2}, ... };
    MSSpectrum(const std::initializer_list<Peak1D>& init);

    /// Copy constructor
    MSSpectrum(const MSSpectrum& source);

    /// Move constructor
    MSSpectrum(MSSpectrum&&) = default;

    /// Destructor
    ~MSSpectrum() = default;

    /// Assignment operator
    MSSpectrum& operator=(const MSSpectrum& source);

    /// Move assignment operator
    MSSpectrum& operator=(MSSpectrum&&) & = default;

    /// Assignment operator
    MSSpectrum& operator=(const SpectrumSettings & source);

    /// Equality operator
    bool operator==(const MSSpectrum& rhs) const;

    /// Equality operator
    bool operator!=(const MSSpectrum& rhs) const
    {
      return !(operator==(rhs));
    }

    // Docu in base class (RangeManager)
    void updateRanges() override;

    ///@name Accessors for meta information
    ///@{
    /// Returns the absolute retention time (in seconds)
    double getRT() const;

    /// Sets the absolute retention time (in seconds)
    void setRT(double rt);

    /**
      @brief Returns the ion mobility drift time (IMTypes::DRIFTTIME_NOT_SET means it is not set)

      @note Drift times may be stored directly as an attribute of the spectrum
      (if they relate to the spectrum as a whole). In case of ion mobility
      spectra, the drift time of the spectrum will always be set here while the
      drift times attribute in the Precursor class may often be unpopulated.
    */
    double getDriftTime() const;

    /**
      @brief Sets the ion mobility drift time
    */
    void setDriftTime(double dt);

    /**
      @brief Returns the ion mobility drift time unit
    */
    DriftTimeUnit getDriftTimeUnit() const;

    /// returns the ion mobility drift time unit as string
    String getDriftTimeUnitAsString() const;

    /**
      @brief Sets the ion mobility drift time unit
    */
    void setDriftTimeUnit(DriftTimeUnit dt);

    /**
      @brief Returns the MS level.

      For survey scans this is 1, for MS/MS scans 2, ...
    */
    UInt getMSLevel() const;

    /// Sets the MS level.
    void setMSLevel(UInt ms_level);

    /// Returns the name
    const String& getName() const;

    /// Sets the name
    void setName(const String& name);

    //@}

    /**
      @name Peak data array methods

      These methods are used to annotate each peak in a spectrum with meta information.
      It is an intermediate way between storing the information in the peak's MetaInfoInterface
      and deriving a new peak type with members for this information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    //@{
    /// Returns a const reference to the float meta data arrays
    const FloatDataArrays& getFloatDataArrays() const;

    /// Returns a mutable reference to the float meta data arrays
    FloatDataArrays& getFloatDataArrays()
    {
      return float_data_arrays_;
    }

    /// Sets the float meta data arrays
    void setFloatDataArrays(const FloatDataArrays& fda);

    /// Returns a const reference to the string meta data arrays
    const StringDataArrays& getStringDataArrays() const;

    /// Returns a mutable reference to the string meta data arrays
    StringDataArrays& getStringDataArrays();

    /// Sets the string meta data arrays
    void setStringDataArrays(const StringDataArrays& sda);

    /// Returns a const reference to the integer meta data arrays
    const IntegerDataArrays& getIntegerDataArrays() const;

    /// Returns a mutable reference to the integer meta data arrays
    IntegerDataArrays& getIntegerDataArrays();

    /// Sets the integer meta data arrays
    void setIntegerDataArrays(const IntegerDataArrays& ida);
    //@}

    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
    */
    void sortByIntensity(bool reverse = false);

    /**
      @brief Lexicographically sorts the peaks by their position.

      The spectrum is sorted with respect to position. Meta data arrays will be sorted accordingly.
    */
    void sortByPosition();

    /**
      @brief Sorts the m/z peaks by their ion mobility value (and the accociated IM data arrays accordingly).

      Requires a binary data array which is a child of 'MS:1002893 ! ion mobility array' (see getIMData())
      
      @throws Exception::MissingInformation if containsIMData() returns false
    */
    void sortByIonMobility();

    /**
      @brief Sort the spectrum, but uses the fact, that certain chunks are presorted
      @param chunks a Chunk is an object that contains the start and end of a sublist of peaks in the spectrum, that is or isn't sorted yet (is_sorted member)
    */
    void sortByPositionPresorted(const std::vector<Chunk>& chunks);

    /// Checks if all peaks are sorted with respect to ascending m/z
    bool isSorted() const;

    /// Checks if m/z peaks are sorted by their associated ion mobility value.
    /// Requires a binary data array which is a child of 'MS:1002893 ! ion mobility array' (see getIMData())
    ///
    /// @throws Exception::MissingInformation if containsIMData() returns false
    bool isSortedByIM() const;

    /// Checks if container is sorted by a certain user-defined property.
    /// You can pass any lambda function with <tt>[](Size index_1, Size index_2) --> bool</tt>
    /// which given two indices into MSSpectrum (either for peaks or data arrays) returns a weak-ordering.
    /// (you need to capture the MSSpectrum in the lambda and operate on it, based on the indices)
    template<class Predicate>
    bool isSorted(const Predicate& lambda) const
    {
      auto value_2_index_wrapper = [this, &lambda](const value_type& value1, const value_type& value2) {
        // translate values into indices (this relies on no copies being made!)
        const Size index1 = (&value1) - (&this->front());
        const Size index2 = (&value2) - (&this->front());
        // just make sure the pointers above are actually pointing to a Peak inside our container
        assert(index1 < this->size()); 
        assert(index2 < this->size());
        return lambda(index1, index2);
      }; 
      return std::is_sorted(this->begin(), this->end(), value_2_index_wrapper);
    }

    /// Sort by a user-defined property
    /// You can pass any @p lambda function with <tt>[](Size index_1, Size index_2) --> bool</tt>
    /// which given two indices into MSSpectrum (either for peaks or data arrays) returns a weak-ordering.
    /// (you need to capture the MSSpectrum in the lambda and operate on it, based on the indices)
    template<class Predicate> 
    void sort(const Predicate& lambda)
    {
      std::vector<Size> indices(this->size());
      std::iota(indices.begin(), indices.end(), 0);
      std::stable_sort(indices.begin(), indices.end(), lambda);
      select(indices);
    }

    //@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific m/z

      @param mz The searched for mass-to-charge ratio searched
      @return Returns the index of the peak.

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the spectrum is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType mz) const;

    /**
      @brief Binary search for the peak nearest to a specific m/z given a +/- tolerance windows in Th

      @param mz The searched for mass-to-charge ratio searched
      @param tolerance The non-negative tolerance applied to both sides of mz

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findNearest(CoordinateType mz, CoordinateType tolerance) const;

    /**
      @brief Search for the peak nearest to a specific m/z given two +/- tolerance windows in Th

      @param mz The searched for mass-to-charge ratio searched
      @param tolerance_left The non-negative tolerance applied left of mz
      @param tolerance_right The non-negative tolerance applied right of mz

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
      @note Search for the left border is done using a binary search followed by a linear scan
    */
    Int findNearest(CoordinateType mz, CoordinateType tolerance_left, CoordinateType tolerance_right) const;

    /**
      @brief Search for the peak with highest intensity among the peaks near to a specific m/z given two +/- tolerance windows in Th

      @param mz The searched for mass-to-charge ratio searched
      @param tolerance_left The non-negative tolerance applied left of mz
      @param tolerance_right The non-negative tolerance applied right of mz
      @return Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findHighestInWindow(CoordinateType mz, CoordinateType tolerance_left, CoordinateType tolerance_right) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator MZBegin(CoordinateType mz);

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator MZEnd(CoordinateType mz);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end);

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator MZBegin(CoordinateType mz) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator MZEnd(CoordinateType mz) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const;

    /**
      @brief Binary search for peak range begin

      Alias for MZBegin()

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator PosBegin(CoordinateType mz);

    /**
      @brief Binary search for peak range begin

      Alias for MZBegin()

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator PosBegin(Iterator begin, CoordinateType mz, Iterator end);

    /**
      @brief Binary search for peak range begin

      Alias for MZBegin()

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator PosBegin(CoordinateType mz) const;

    /**
      @brief Binary search for peak range begin

      Alias for MZBegin()

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator PosBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MZEnd()

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator PosEnd(CoordinateType mz);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MZEnd()

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator PosEnd(Iterator begin, CoordinateType mz, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MZEnd()

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator PosEnd(CoordinateType mz) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for MZEnd()

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator PosEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const;

    /// do the names of internal float metadata arrays contain any hint of ion mobility data, i.e. they are a child of 'MS:1002893 ! ion mobility array'?
    /// (for spectra which represent an IM-frame)
    bool containsIMData() const;

    /**
      @brief Get the Ion mobility data array's @p index and its associated @p unit

      This only works for spectra which represent an IM-frame, i.e. they have a float metadata array which is a child of 'MS:1002893 ! ion mobility array'?
      Check this first by using `containsIMData()`.

      @throws Exception::MissingInformation if IM data is not present
    */
    std::pair<Size, DriftTimeUnit> getIMData() const;
    

    /**
      @brief Get the spectrum's ion mobility data (if exists) and its associated unit as a pair of {unit, data}
      This only works for spectra which represent an IM-frame, i.e. they have a float metadata array which is a child of 'MS:1002893 ! ion mobility array'.
      If this is not present, this returns {DriftTimeUnit::NONE, {}}
    */
    std::pair<DriftTimeUnit, std::vector<float>> maybeGetIMData() const;
    
    //@}


    /**
      @brief Clears all data and meta data

      Will delete (clear) all peaks contained in the spectrum as well as any
      associated data arrays (FloatDataArrays, IntegerDataArrays,
      StringDataArrays) by default. If @em clear_meta_data is @em true, then
      also all meta data (such as RT, drift time, ms level etc) will be
      deleted.

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    /*
      @brief Select a (subset of) spectrum and its data_arrays, only retaining the indices given in @p indices

      @param indices Vector of indices to keep
      @return Reference to this MSSpectrum

    */
    MSSpectrum& select(const std::vector<Size>& indices);


    /**
      @brief Determine if spectrum is profile or centroided using up to three layers of information.

      First, the SpectrumSettings are inquired and the type is returned unless it is unknown.
      Second, all data processing entries are searched for a centroiding step.
      If that is unsuccessful as well and @p query_data is true, the data is fed into PeakTypeEstimator().

      @param [query_data] If SpectrumSettings and DataProcessing information are not sufficient, should the data be queried? (potentially expensive)
      @return The spectrum type (centroided, profile or unknown)
    */
    SpectrumSettings::SpectrumType getType(const bool query_data) const;
    using SpectrumSettings::getType; // expose base class function

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the spectrum is unsorted.
    ConstIterator getBasePeak() const;

    /// return the peak with the highest intensity. If the peak is not unique, the first peak in the container is returned.
    /// The function works correctly, even if the spectrum is unsorted.
    Iterator getBasePeak();

    /// compute the total ion count (sum of all peak intensities)
    PeakType::IntensityType calculateTIC() const;

protected:
    /// Retention time
    double retention_time_ = -1;

    /// Drift time
    double drift_time_ = -1;

    /// Drift time unit
    DriftTimeUnit drift_time_unit_ = DriftTimeUnit::NONE;

    /// MS level
    UInt ms_level_ = 1;

    /// Name
    String name_;

    /// Float data arrays
    FloatDataArrays float_data_arrays_;

    /// String data arrays
    StringDataArrays string_data_arrays_;

    /// Integer data arrays
    IntegerDataArrays integer_data_arrays_;
  };

  inline std::ostream& operator<<(std::ostream& os, const MSSpectrum& spec)
  {
    os << "-- MSSPECTRUM BEGIN --" << std::endl;

    // spectrum settings
    os << static_cast<const SpectrumSettings&>(spec);

    // peaklist
    for (MSSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
    {
      os << *it << std::endl;
    }

    os << "-- MSSPECTRUM END --" << std::endl;
    return os;
  }

} // namespace OpenMS
