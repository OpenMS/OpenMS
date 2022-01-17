// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <cstddef> // for ptrdiff_t
#include <string>
#include <vector>


namespace OpenMS
{
  /**
    @brief Class to hold strings, numeric values, vectors of strings and vectors of numeric values using the stl types.

    - To choose one of these types, just use the appropriate constructor.
    - Automatic conversion is supported and throws Exceptions in case of invalid conversions.
    - An empty object is created with the default constructor.

    @ingroup Datastructures
  */
  class OPENMS_DLLAPI ParamValue
  {

public:

    /// Empty data value for comparisons
    static const ParamValue EMPTY;

    /// Supported types for ParamValue
    enum ValueType : unsigned char
    {
      STRING_VALUE, ///< string value
      INT_VALUE, ///< integer value
      DOUBLE_VALUE, ///< double value
      STRING_LIST, ///< string vector
      INT_LIST, ///< integer vector
      DOUBLE_LIST, ///< double vector
      EMPTY_VALUE ///< empty value
    };

    /// @name Constructors and destructors
    //@{
    /// Default constructor
    ParamValue();
    /// Copy constructor
    ParamValue(const ParamValue&);
    /// Move constructor
    ParamValue(ParamValue&&) noexcept;
    /// specific constructor for char* (converted to string)
    ParamValue(const char*);
    /// specific constructor for std::string values
    ParamValue(const std::string&);
    /// specific constructor for string vectors
    ParamValue(const std::vector<std::string>&);
    /// specific constructor for integer vectors
    ParamValue(const std::vector<int>&);
    /// specific constructor for double vectors
    ParamValue(const std::vector<double>&);
    /// specific constructor for long double values (note: the implementation uses double)
    ParamValue(long double);
    /// specific constructor for double values (note: the implementation uses double)
    ParamValue(double);
    /// specific constructor for float values (note: the implementation uses double)
    ParamValue(float);
    /// specific constructor for short int values (note: the implementation uses ptrdiff_t)
    ParamValue(short int);
    /// specific constructor for unsigned short int values (note: the implementation uses ptrdiff_t)
    ParamValue(unsigned short int);
    /// specific constructor for int values (note: the implementation uses ptrdiff_t)
    ParamValue(int);
    /// specific constructor for unsigned int values (note: the implementation uses ptrdiff_t)
    ParamValue(unsigned);
    /// specific constructor for long int values (note: the implementation uses ptrdiff_t)
    ParamValue(long int);
    /// specific constructor for unsigned long int values (note: the implementation uses ptrdiff_t)
    ParamValue(unsigned long);
    /// specific constructor for long long int values (note: the implementation uses ptrdiff_t)
    ParamValue(long long);
    /// specific constructor for unsigned long long int values (note: the implementation uses ptrdiff_t)
    ParamValue(unsigned long long);
    /// Destructor
    ~ParamValue();
    //@}

    ///@name Cast operators
    ///These methods are used when the DataType is known.
    ///If they are applied to a ParamValue with the wrong DataType, an exception (Exception::ConversionError) is thrown. In particular, none of these operators will work for an empty ParamValue (DataType EMPTY_VALUE) - except toChar(), which will return 0.
    //@{

    /**
      @brief conversion operator to string

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator std::string() const;

    /**
      @brief conversion operator to string vector

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator std::vector<std::string>() const;

    /**
      @brief conversion operator to integer vector

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator std::vector<int>() const;

    /**
      @brief conversion operator to double vector

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator std::vector<double>() const;

    /**
      @brief conversion operator to long double

      Note: The implementation uses typedef double (as opposed to float, double, long double.)

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator long double() const;

    /**
      @brief conversion operator to double

      Note: The implementation uses typedef double (as opposed to float, double, long double.)

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator double() const;

    /**
      @brief conversion operator to float

      Note: The implementation uses typedef double (as opposed to float, double, long double.)

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator float() const;

    /**
      @brief conversion operator to short int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator short int() const;

    /**
      @brief conversion operator to unsigned short int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned short int() const;

    /**
      @brief conversion operator to int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */

    operator int() const;

    /**
      @brief conversion operator to unsigned int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned int() const;

    /**
      @brief conversion operator to long int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator long int() const;

    /**
      @brief conversion operator to unsigned long int

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned long int() const;

    /**
      @brief conversion operator to long long

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator long long() const;

    /**
      @brief conversion operator to unsigned long long

      Note: The implementation uses typedef ptrdiff_t.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned long long() const;

    /**
      @brief Conversion to bool

      Converts the strings 'true' and 'false' to a bool.

      @exception Exception::ConversionError is thrown for non-string parameters and string parameters with values other than 'true' and 'false'.
    */
    bool toBool() const;

    /**
      @brief Convert ParamValues to char*

      If the ParamValue contains a string, a pointer to it's char* is returned.
      If the ParamValue is empty, nullptr is returned.
    */
    const char* toChar() const;

    /**
     * @brief Convert ParamValue to string
     *
     * @exception Exception::ConversionError is thrown for ParamValue::EMPTY and
     */

    std::string toString(bool full_precision = true) const;

    /**
      @brief Explicitly convert ParamValue to string vector

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    std::vector<std::string> toStringVector() const;

    /**
      @brief Explicitly convert ParamValue to IntList

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    std::vector<int> toIntVector() const;

    /**
      @brief Explicitly convert ParamValue to DoubleList

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    std::vector<double> toDoubleVector() const;
    //@}

    ///@name Assignment operators
    ///These methods are used to assign supported types directly to a ParamValue object.
    //@{
    /// Assignment operator
    ParamValue& operator=(const ParamValue&);
    /// Move assignment operator
    ParamValue& operator=(ParamValue&&) noexcept;
    /// specific assignment for char* (converted to string)
    ParamValue& operator=(const char*);
    /// specific assignment for std::string values
    ParamValue& operator=(const std::string&);
    /// specific assignment for string vectors
    ParamValue& operator=(const std::vector<std::string>&);
    /// specific assignment for integer vectors
    ParamValue& operator=(const std::vector<int>&);
    /// specific assignment for double vectors
    ParamValue& operator=(const std::vector<double>&);
    /// specific assignment for long double values (note: the implementation uses double)
    ParamValue& operator=(const long double);
    /// specific assignment for double values (note: the implementation uses double)
    ParamValue& operator=(const double);
    /// specific assignment for float values (note: the implementation uses double)
    ParamValue& operator=(const float);
    /// specific assignment for short int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const short int);
    /// specific assignment for unsigned short int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const unsigned short int);
    /// specific assignment for int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const int);
    /// specific assignment for unsigned int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const unsigned);
    /// specific assignment for long int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const long int);
    /// specific assignment for unsigned long int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const unsigned long);
    /// specific assignment for long long int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const long long);
    /// specific assignment for unsigned long long int values (note: the implementation uses ptrdiff_t)
    ParamValue& operator=(const unsigned long long);
    //@}

    /// returns the type of value stored
    inline ValueType valueType() const
    {
      return value_type_;
    }

    /**
       @brief Test if the value is empty

       @note A ParamValue containing an empty string ("") does not count as empty!
    */
    inline bool isEmpty() const
    {
      return value_type_ == EMPTY_VALUE;
    }

    /// output stream operator
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream&, const ParamValue&);

    /// Equality comparator
    friend OPENMS_DLLAPI bool operator==(const ParamValue&, const ParamValue&);

    /// Smaller than comparator (for vectors we use the size)
    friend OPENMS_DLLAPI bool operator<(const ParamValue&, const ParamValue&);

    /// Greater than comparator (for vectors we use the size)
    friend OPENMS_DLLAPI bool operator>(const ParamValue&, const ParamValue&);

    /// Equality comparator
    friend OPENMS_DLLAPI bool operator!=(const ParamValue&, const ParamValue&);

protected:

    /// Type of the currently stored value
    ValueType value_type_;

    /// Space to store the data
    union
    {
      ptrdiff_t ssize_;
      double dou_;
      std::string* str_;
      std::vector<std::string>* str_list_;
      std::vector<int>* int_list_;
      std::vector<double>* dou_list_;
    } data_;

private:

    /// Clears the current state of the ParamValue and release every used memory.
    void clear_() noexcept;

    /// Convert a double to std::string
    /// with full precision 15 decimal places are given, otherwise 3
    /// numbers above 10000 or below 0.0001 are given in scientific notation (i.e. 1.0e04)
    static std::string doubleToString(double value, bool full_precision = true);

  };
}

