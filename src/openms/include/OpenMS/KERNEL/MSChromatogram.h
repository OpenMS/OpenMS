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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSCHROMATOGRAM_H
#define OPENMS_KERNEL_MSCHROMATOGRAM_H

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/METADATA/ChromatogramSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/METADATA/DataArrays.h>

namespace OpenMS
{
  class ChromatogramPeak;

  /**
    @brief The representation of a chromatogram.
    @ingroup Kernel
  */

  class OPENMS_DLLAPI MSChromatogram :
    private std::vector<ChromatogramPeak>,
    public RangeManager<1>,
    public ChromatogramSettings
  {

public:

    /// Comparator for the retention time.
    struct OPENMS_DLLAPI MZLess :
      public std::binary_function<MSChromatogram, MSChromatogram, bool>
    {
      bool operator()(const MSChromatogram& a, const MSChromatogram& b) const;
    };

    ///@name Base type definitions
    ///@{
    /// Peak type
    typedef ChromatogramPeak PeakType;
    /// Coordinate (RT) type
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Chromatogram base type
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

    ///@name Export methods from std::vector
    //@{
    using ContainerType::operator[];
    using ContainerType::begin;
    using ContainerType::rbegin;
    using ContainerType::end;
    using ContainerType::rend;
    using ContainerType::resize;
    using ContainerType::size;
    using ContainerType::push_back;
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
    MSChromatogram() :
      ContainerType(),
      RangeManager<1>(),
      ChromatogramSettings(),
      name_(),
      float_data_arrays_(),
      string_data_arrays_(),
      integer_data_arrays_()
    {}

    /// Copy constructor
    MSChromatogram(const MSChromatogram& source) :
      ContainerType(source),
      RangeManager<1>(source),
      ChromatogramSettings(source),
      name_(source.name_),
      float_data_arrays_(source.float_data_arrays_),
      string_data_arrays_(source.string_data_arrays_),
      integer_data_arrays_(source.integer_data_arrays_)
    {}

    /// Destructor
    ~MSChromatogram() override
    {}

    /// Assignment operator
    MSChromatogram& operator=(const MSChromatogram& source);

    /// Equality operator
    bool operator==(const MSChromatogram& rhs) const;

    /// Equality operator
    bool operator!=(const MSChromatogram& rhs) const
    {
      return !(operator==(rhs));
    }

    // Docu in base class (RangeManager)
    void updateRanges() override
    {
      this->clearRanges();
      updateRanges_(ContainerType::begin(), ContainerType::end());
    }

    ///@name Accessors for meta information
    ///@{
    /// Returns the name
    const String& getName() const;

    /// Sets the name
    void setName(const String& name);

    ///@}

    /// returns the mz of the product entry, makes sense especially for MRM scans
    double getMZ() const;

    /**
      @name Peak data array methods

      These methods are used to annotate each peak in a chromatogram with meta information.
      It is an intermediate way between storing the information in the peak's MetaInfoInterface
      and deriving a new peak type with members for this information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    ///@{
    /// Returns a const reference to the float meta data arrays
    const FloatDataArrays& getFloatDataArrays() const;

    /// Returns a mutable reference to the float meta data arrays
    FloatDataArrays& getFloatDataArrays();

    /// Sets the float meta data arrays
    void setFloatDataArrays(const FloatDataArrays& fda)
    {
      float_data_arrays_ = fda;
    }

    /// Returns a const reference to the string meta data arrays
    const StringDataArrays& getStringDataArrays() const;

    /// Returns a mutable reference to the string meta data arrays
    StringDataArrays& getStringDataArrays();

    /// Sets the string meta data arrays
    void setStringDataArrays(const StringDataArrays& sda)
    {
      string_data_arrays_ = sda;
    }

    /// Returns a const reference to the integer meta data arrays
    const IntegerDataArrays& getIntegerDataArrays() const;

    /// Returns a mutable reference to the integer meta data arrays
    IntegerDataArrays& getIntegerDataArrays();

    /// Sets the integer meta data arrays
    void setIntegerDataArrays(const IntegerDataArrays& ida)
    {
      integer_data_arrays_ = ida;
    }

    ///@}

    ///@name Sorting peaks
    ///@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
    */
    void sortByIntensity(bool reverse = false);

    /**
      @brief Lexicographically sorts the peaks by their position.

      The chromatogram is sorted with respect to position. Meta data arrays will be sorted
      accordingly.
    */
    void sortByPosition();

    ///Checks if all peaks are sorted with respect to ascending RT
    bool isSorted() const;

    ///@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific RT

      @param rt The searched for mass-to-charge ratio searched
      @return Returns the index of the peak.

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the chromatogram is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType rt) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to retention time! Otherwise the
      result is undefined.
    */
    Iterator RTBegin(CoordinateType rt);

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    Iterator RTBegin(Iterator begin, CoordinateType rt, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator RTEnd(CoordinateType rt);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator RTEnd(Iterator begin, CoordinateType rt, Iterator end);

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator RTBegin(CoordinateType rt) const;

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator RTBegin(ConstIterator begin, CoordinateType rt, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator RTEnd(CoordinateType rt) const;

    ConstIterator MZEnd(CoordinateType rt) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator RTEnd(ConstIterator begin, CoordinateType rt, ConstIterator end) const;

    /**
      @brief Binary search for peak range begin

      Alias for RTBegin()

      @note Make sure the chromatogram is sorted with respect to retention time! Otherwise the
      result is undefined.
    */
    Iterator PosBegin(CoordinateType rt);

    /**
      @brief Binary search for peak range begin

      Alias for RTBegin()

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    Iterator PosBegin(Iterator begin, CoordinateType rt, Iterator end);

    /**
      @brief Binary search for peak range begin

      Alias for RTBegin()

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator PosBegin(CoordinateType rt) const;

    /**
      @brief Binary search for peak range begin

      Alias for RTBegin()

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator PosBegin(ConstIterator begin, CoordinateType rt, ConstIterator end) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for RTEnd()

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator PosEnd(CoordinateType rt);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for RTEnd()

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator PosEnd(Iterator begin, CoordinateType rt, Iterator end);

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for RTEnd()

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator PosEnd(CoordinateType rt) const;

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      Alias for RTEnd()

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator PosEnd(ConstIterator begin, CoordinateType rt, ConstIterator end) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    ///@}

protected:

    /// Name
    String name_;

    /// Float data arrays
    FloatDataArrays float_data_arrays_;

    /// String data arrays
    StringDataArrays string_data_arrays_;

    /// Integer data arrays
    IntegerDataArrays integer_data_arrays_;
  };

  /// Print the contents to a stream.
  std::ostream& operator<<(std::ostream& os, const MSChromatogram& chrom);

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSCHROMATOGRAM_H
