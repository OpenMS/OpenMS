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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

namespace OpenMS
{
  class Peak1D;
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
  class OPENMS_DLLAPI MSSpectrum :
    private std::vector<Peak1D>,
    public RangeManager<1>,
    public SpectrumSettings
  {
public:

    /// Comparator for the retention time.
    struct OPENMS_DLLAPI RTLess : public std::binary_function<MSSpectrum, MSSpectrum, bool>
    {
      bool operator()(const MSSpectrum& a, const MSSpectrum& b) const;
    };

    ///@name Base type definitions
    //@{
    /// Peak type
    typedef OpenMS::Peak1D PeakType;
    /// Coordinate (m/z) type
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Spectrum base type
    typedef std::vector<PeakType> ContainerType;
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

    /// Copy constructor
    MSSpectrum(const MSSpectrum& source);

    /// Destructor
    ~MSSpectrum() override
    {}

    /// Assignment operator
    MSSpectrum& operator=(const MSSpectrum& source);

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
      @brief Returns the ion mobility drift time in milliseconds (-1 means it is not set)

      @note Drift times may be stored directly as an attribute of the spectrum
      (if they relate to the spectrum as a whole). In case of ion mobility
      spectra, the drift time of the spectrum will always be set here while the
      drift times attribute in the Precursor class may often be unpopulated.
    */
    double getDriftTime() const;

    /**
      @brief Returns the ion mobility drift time in milliseconds
    */
    void setDriftTime(double dt);

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

    /// Returns a mutable reference to the first integer meta data array with the given name
    inline IntegerDataArray& getIntegerDataArrayByName(String name)
    {
      return *std::find_if(integer_data_arrays_.begin(), integer_data_arrays_.end(), 
        [&name](const IntegerDataArray& da) { return da.getName() == name; } );
    }

    /// Returns a mutable reference to the first string meta data array with the given name
    inline StringDataArray& getStringDataArrayByName(String name)
    {
      return *std::find_if(string_data_arrays_.begin(), string_data_arrays_.end(), 
        [&name](const StringDataArray& da) { return da.getName() == name; } );
    }

    /// Returns a mutable reference to the first float meta data array with the given name
    inline FloatDataArray& getFloatDataArrayByName(String name)
    {
      return *std::find_if(float_data_arrays_.begin(), float_data_arrays_.end(), 
        [&name](const FloatDataArray& da) { return da.getName() == name; } );
    }

    /// Returns a const reference to the first integer meta data array with the given name
    inline const IntegerDataArray& getIntegerDataArrayByName(String name) const
    {
      return *std::find_if(integer_data_arrays_.begin(), integer_data_arrays_.end(), 
        [&name](const IntegerDataArray& da) { return da.getName() == name; } );
    }

    /// Returns a const reference to the first string meta data array with the given name
    inline const StringDataArray& getStringDataArrayByName(String name) const
    {
      return *std::find_if(string_data_arrays_.begin(), string_data_arrays_.end(), 
        [&name](const StringDataArray& da) { return da.getName() == name; } );
    }

    /// Returns a const reference to the first float meta data array with the given name
    inline const FloatDataArray& getFloatDataArrayByName(String name) const
    {
      return *std::find_if(float_data_arrays_.begin(), float_data_arrays_.end(), 
        [&name](const FloatDataArray& da) { return da.getName() == name; } );
    }

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

    /// Checks if all peaks are sorted with respect to ascending m/z
    bool isSorted() const;

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

    //@}


    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    /*
      @brief Select a (subset of) spectrum and its data_arrays, only retaining the indices given in @p indices

      @param indices Vector of indices to keep
      @return Reference to this MSSpectrum

    */
    MSSpectrum& select(const std::vector<Size>& indices);

protected:

    /// Retention time
    double retention_time_;

    /// Drift time
    double drift_time_;

    /// MS level
    UInt ms_level_;

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

#endif // OPENMS_KERNEL_MSSPECTRUM_H
