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

#ifndef OPENMS_DATASTRUCTURES_DATAVALUE_H
#define OPENMS_DATASTRUCTURES_DATAVALUE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

class QString;

namespace OpenMS
{
  /**
    @brief Class to hold strings, numeric values, lists of strings and lists of numeric values.

    - To choose one of these types, just use the appropriate constructor.
    - Automatic conversion is supported and throws Exceptions in case of invalid conversions.
    - An empty object is created with the default constructor.

    @ingroup Datastructures
  */
  class OPENMS_DLLAPI DataValue
  {

public:

    /// Empty data value for comparisons
    static const DataValue EMPTY;

    /// Supported types for DataValue
    enum DataType
    {
      STRING_VALUE, ///< string value
      INT_VALUE, ///< integer value
      DOUBLE_VALUE, ///< double value
      STRING_LIST, ///< string list
      INT_LIST, ///< integer list
      DOUBLE_LIST, ///< double list
      EMPTY_VALUE ///< empty value
    };

    /// @name Constructors and destructors
    //@{
    /// default constructor
    DataValue();
    /// specific constructor for char* (converted to string)
    DataValue(const char*);
    /// specific constructor for std::string values
    DataValue(const std::string&);
    /// specific constructor for string values
    DataValue(const String&);
    /// specific constructor for QString values
    DataValue(const QString&);
    /// specific constructor for string lists
    DataValue(const StringList&);
    /// specific constructor for integer lists
    DataValue(const IntList&);
    /// specific constructor for double lists
    DataValue(const DoubleList&);
    /// specific constructor for long double values (note: the implementation uses double)
    DataValue(long double);
    /// specific constructor for double values (note: the implementation uses double)
    DataValue(double);
    /// specific constructor for float values (note: the implementation uses double)
    DataValue(float);
    /// specific constructor for short int values (note: the implementation uses SignedSize)
    DataValue(short int);
    /// specific constructor for unsigned short int values (note: the implementation uses SignedSize)
    DataValue(unsigned short int);
    /// specific constructor for int values (note: the implementation uses SignedSize)
    DataValue(int);
    /// specific constructor for unsigned int values (note: the implementation uses SignedSize)
    DataValue(unsigned);
    /// specific constructor for long int values (note: the implementation uses SignedSize)
    DataValue(long int);
    /// specific constructor for unsigned long int values (note: the implementation uses SignedSize)
    DataValue(unsigned long);
    /// specific constructor for long long int values (note: the implementation uses SignedSize)
    DataValue(long long);
    /// specific constructor for unsigned long long int values (note: the implementation uses SignedSize)
    DataValue(unsigned long long);
    /// copy constructor
    DataValue(const DataValue&);
    /// destructor
    virtual ~DataValue();
    //@}

    ///@name cast operators
    ///These methods are used when the DataType is known.
    ///If they are applied to a DataValue with the wrong DataType, an exception is thrown.
    //@{

    /**
      @brief conversion operator to string

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator std::string() const;

    /**
      @brief conversion operator to string list

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator StringList() const;

    /**
      @brief conversion operator to integer list

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator IntList() const;

    /**
      @brief conversion operator to double list

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator DoubleList() const;

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

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator short int() const;

    /**
      @brief conversion operator to unsigned short int

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned short int() const;

    /**
      @brief conversion operator to int

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */

    operator int() const;

    /**
      @brief conversion operator to unsigned int

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned int() const;

    /**
      @brief conversion operator to long int

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator long int() const;

    /**
      @brief conversion operator to unsigned long int

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned long int() const;

    /**
      @brief conversion operator to long long

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator long long() const;

    /**
      @brief conversion operator to unsigned long long

      Note: The implementation uses typedef SignedSize.

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    operator unsigned long long() const;

    /**
      @brief Convert DataValues to char*

      If the DataValue contains a string, a pointer to it's char* is returned.
      If the DataValue is empty, NULL is returned.
    */
    const char* toChar() const;

    /**
      @brief Explicitly convert DataValue to StringList

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    StringList toStringList() const;

    /**
      @brief Explicitly convert DataValue to IntList

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    IntList toIntList() const;

    /**
      @brief Explicitly convert DataValue to DoubleList

      @exception Exception::ConversionError is thrown if a cast from the the wrong type is requested
    */
    DoubleList toDoubleList() const;
    //@}

    ///@name assignment/conversion operators
    ///These methods are used to assign supported types to DataType.
    //@{

    /// specific assignment for char* (converted to string)
    DataValue& operator=(const char*);
    /// specific assignment for std::string values
    DataValue& operator=(const std::string&);
    /// specific assignment for string values
    DataValue& operator=(const String&);
    /// specific assignment for QString values
    DataValue& operator=(const QString&);
    /// specific assignment for string lists
    DataValue& operator=(const StringList&);
    /// specific assignment for integer lists
    DataValue& operator=(const IntList&);
    /// specific assignment for double lists
    DataValue& operator=(const DoubleList&);
    /// specific assignment for long double values (note: the implementation uses double)
    DataValue& operator=(const long double);
    /// specific assignment for double values (note: the implementation uses double)
    DataValue& operator=(const double);
    /// specific assignment for float values (note: the implementation uses double)
    DataValue& operator=(const float);
    /// specific assignment for short int values (note: the implementation uses SignedSize)
    DataValue& operator=(const short int);
    /// specific assignment for unsigned short int values (note: the implementation uses SignedSize)
    DataValue& operator=(const unsigned short int);
    /// specific assignment for int values (note: the implementation uses SignedSize)
    DataValue& operator=(const int);
    /// specific assignment for unsigned int values (note: the implementation uses SignedSize)
    DataValue& operator=(const unsigned);
    /// specific assignment for long int values (note: the implementation uses SignedSize)
    DataValue& operator=(const long int);
    /// specific assignment for unsigned long int values (note: the implementation uses SignedSize)
    DataValue& operator=(const unsigned long);
    /// specific assignment for long long int values (note: the implementation uses SignedSize)
    DataValue& operator=(const long long);
    /// specific assignment for unsigned long long int values (note: the implementation uses SignedSize)
    DataValue& operator=(const unsigned long long);

    //@}

    ///@name conversion operators
    ///These methods can be used independent of the DataType. If you already know the DataType, you should use a cast operator!
    /// <BR>For conversion of string DataValues to numeric types, first use toString() and then the conversion methods of String.
    //@{

    ///Conversion to String.
    String toString() const;

    ///Conversion to QString
    QString toQString() const;

    /**
      @brief Conversion to bool

      Converts the strings 'true' and 'false' to a bool.

      @exception Exception::ConversionError is thrown for non-string parameters and string parameters with values other than 'true' and 'false'.
    */
    bool toBool() const;
    //@}

    /// returns the type of value stored
    inline DataType valueType() const
    {
      return value_type_;
    }

    /// assignment operator
    DataValue& operator=(const DataValue&);

    /// test if the value is empty
    inline bool isEmpty() const
    {
      return value_type_ == EMPTY_VALUE;
    }

    ///@name Methods to handle units
    ///These methods are used when the DataValue has an associated unit.
    //@{

    /// Check if the value has a unit
    inline bool hasUnit() const
    {
      return unit_ != "";
    }

    /// Return the unit associated to this DataValue.
    const String& getUnit() const;

    /// Sets the unit to the given String.
    void setUnit(const String& unit);

    //@}

    /// output stream operator
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream&, const DataValue&);

    /// Equality comparator
    friend OPENMS_DLLAPI bool operator==(const DataValue&, const DataValue&);

    /// Smaller than comparator (for lists we use the size)
    friend OPENMS_DLLAPI bool operator<(const DataValue&, const DataValue&);

    /// Greater than comparator (for lists we use the size)
    friend OPENMS_DLLAPI bool operator>(const DataValue&, const DataValue&);

    /// Equality comparator
    friend OPENMS_DLLAPI bool operator!=(const DataValue&, const DataValue&);

protected:

    /// Type of the currently stored value
    DataType value_type_;

    /// Space to store the data
    union
    {
      SignedSize ssize_;
      double dou_;
      String* str_;
      StringList* str_list_;
      IntList* int_list_;
      DoubleList* dou_list_;
    } data_;

private:
    /// The unit of the data value (if it has one), otherwise empty string.
    String unit_;

    /// Clears the current state of the DataValue and release every used memory.
    void clear_();
  };
}

#endif // OPENMS_DATASTRUCTURES_DATAVALUE_H
