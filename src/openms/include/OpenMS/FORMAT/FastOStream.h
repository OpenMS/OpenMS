// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: $
// $Authors: Tom Lukas Lankenau, Anton Kriese $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <vector>
#include <ostream>
#include <cstring>

namespace OpenMS
{
  /**
      @brief Stream wrapper class that writes without using sentries.

      Usage: Construct with an existing std::ostream (or ostream-derived class). That ostream is then referenced by the member os_.

      The operator << () and the function write() make use of the stream buffer-function sputn(), which writes to the stream without
      constructing a sentry on it. This method is unsave if several threads access the stream, but in most cases, only one thread is
      actually accessing the stream. Avoiding the sentry construction and destruction saves runtime and makes the process of writing
      files faster.

      The optimization is being used for: OpenMS::String, std::string, char*, arithmetic types, std::vector, OpenMS::DataValue.

      All other types are written to the ostream via the ostream::operator <<().

      For arithmetic types, String::append(), which calls boost-functions, is used. This gives the best relative speedup amongst all listed types.

      @ingroup FileIO
  */
  class FastOStream
  {
public:
    // Rule of six
    /// Default constructor
    FastOStream() = delete;
    FastOStream(const FastOStream& rhs) = delete;
    FastOStream& operator=(const FastOStream& rhs) = delete;
    FastOStream(FastOStream&& rhs) = delete;
    FastOStream& operator=(FastOStream&& rhs) = delete;
    /// Destructor
    ~FastOStream() = default;

    /**
         @brief Constructor

         @param os The ostream wrapped by the class; Assigned to the member ostream-reference
    */
    FastOStream(std::ostream& os) : os_(os) {};

    /**
         @brief overloads for the operator << ()

         @param s Any kind of string (std::string, OpenMS::String, QString, char*)

         @returns reference to this
    */
    inline FastOStream& operator << (const OpenMS::String& s)
    {
      os_.rdbuf()->sputn(s.c_str(), s.size());
      return *this;
    }

    inline FastOStream& operator << (const std::string& s)
    {
      os_.rdbuf()->sputn(s.c_str(), s.size());
      return *this;
    }

    inline FastOStream& operator << (const char* const s)
    {
      os_.rdbuf()->sputn(s, strlen(s));
      return *this;
    }

    /**
         @brief template overload for operator << ()

         @param v std::vector of type T
         @returns reference to this
    */
    template <typename T>
    inline FastOStream& operator << (const std::vector<T>& v)
    {
      *this << "[";

      if (!v.empty())
      {
        auto end = v.end() - 1;
        for (auto it = v.begin(); it < end; ++it)
        { // convert T to String using FastOStream::operator<< (very fast!)
          *this << *it << ", ";
        }
        *this << v.back();
      }

      *this << "]";
      return *this;
    }

    /**
         @brief template overload for operator << ()

         Arithmetic types are converted to OpenMS::String using boost-functions (in String.append())
         All other types are put into the ostream via the default ostream::operator << ()

         @param s Any other type
         @returns reference to this
    */
    // if constexpr() IS ONLY ALLOWED FROM c++17!
    //template <typename T>
    //inline FastOStream& operator << (const T& s)
    //{
    //  if constexpr (std::is_arithmetic<typename std::decay<T>::type>::value)
    //  {
    //    buffer_.clear();
    //    buffer_ += s; // use internal buffer, instead of String(s), which would allocate
    //    os_.rdbuf()->sputn(buffer_.c_str(), buffer_.size());
    //  }
    //  else
    //  {
    //    os_ << s;
    //  }
    //  return *this;
    //}

    // if input type is arithmetic
    template <typename T>
    inline typename std::enable_if<std::is_arithmetic<typename std::decay<T>::type>::value, FastOStream&>::type
    operator << (const T& value)
    {
      buffer_.clear();
      buffer_ += value; // use internal buffer, instead of String(s), which would allocate
      os_.rdbuf()->sputn(buffer_.c_str(), buffer_.size());
      return *this;
    }

    // if input type is NOT arithmetic
    template <typename T>
    inline typename std::enable_if<!std::is_arithmetic<typename std::decay<T>::type>::value, FastOStream&>::type
    operator << (const T& value)
    {
      os_ << value;
      return *this;
    }

    /**
         @brief overloads for write()

         @param s Any kind of string (std::string, OpenMS::String, QString, char*)
         @param len Length to be written into the stream.
    */
    inline void write(const OpenMS::String& s, uint64_t len)
    {
      write(s.c_str(), len);
    }

    inline void write(const std::string& s, uint64_t len)
    {
      write(s.c_str(), len);
    }

    inline void write(const char* const s, uint64_t len)
    {
      uint64_t written = os_.rdbuf()->sputn(s, len);
      if (written != len) os_.setstate(std::ios_base::badbit);
    }

    /**
         @brief Get-Function to access the wrapped ostream

         @returns reference to the ostream
    */
    std::ostream& getStream()
    {
      return os_;
    }

private:
    std::ostream& os_; ///< Reference to the ostream that is written to
    String buffer_; ///< String, used to convert arithmetic values into String
  };


  /**
       @brief Free function overloading FastOStream::operator<<(const DataValue&)

       @param os The FastOStream to write the DataValue to
       @param p The DataValue to write to the FOS

       @returns reference to FastOStream
  */
  inline FastOStream& operator << (FastOStream& os, const DataValue& p)
  {
    /// for doubles or lists of doubles, you get full precision. Use DataValue::toString(false) if you only need low precision

    switch (p.value_type_)
    {
      case DataValue::STRING_VALUE: os << *(p.data_.str_); break;
      case DataValue::STRING_LIST: os << *(p.data_.str_list_); break;
      case DataValue::INT_LIST: os << *(p.data_.int_list_); break;
      case DataValue::DOUBLE_LIST: os << *(p.data_.dou_list_); break;
      case DataValue::INT_VALUE: os << p.data_.ssize_; break; // using our String conversion (faster than std::ofstream)
      case DataValue::DOUBLE_VALUE: os << p.data_.dou_; break; // using our String conversion (faster than std::ofstream)
      case DataValue::EMPTY_VALUE: break;
    }
    return os;
  }
}
