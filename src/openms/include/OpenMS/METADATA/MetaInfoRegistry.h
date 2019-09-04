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
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <map>
#include <string>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <unordered_map>

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( push )
#pragma warning( disable : 4251 )     // disable MSVC dll-interface warning
#endif

namespace OpenMS
{

  /**
      @brief Registry which assigns unique integer indices to strings.

      When registering a new name an index >= 1024 is assigned.
      Indices from 1 to 1023 are reserved for fast access and will never change:<BR>
      1 - isotopic_range<BR>
      2 - cluster_id<BR>
      3 - label<BR>
      4 - icon<BR>
      5 - color<BR>
      6 - RT<BR>
      7 - MZ<BR>
      8 - predicted_RT<BR>
      9 - predicted_RT_p_value<BR>
      10 - spectrum_reference<BR>
      11 - ID<BR>
      12 - low_quality<BR>
      13 - charge<BR>

      @ingroup Metadata
  */
  class OPENMS_DLLAPI MetaInfoRegistry
  {
public:
    /// Default constructor
    MetaInfoRegistry();

    /// Copy constructor
    MetaInfoRegistry(const MetaInfoRegistry& rhs);

    /// Destructor
    ~MetaInfoRegistry();

    /// Assignment operator
    MetaInfoRegistry& operator=(const MetaInfoRegistry& rhs);

    /**
        Registers a string, stores its description and unit, and returns the corresponding index.
        If the string is already registered, it returns the index of the string.
    */
    UInt registerName(const String& name, const String& description = "", const String& unit = "");

    /**
      @brief Sets the description (String), corresponding to an index

      @exception Exception::InvalidValue is thrown for unregistered indices
    */
    void setDescription(UInt index, const String& description);

    /**
      @brief Sets the description (String), corresponding to a name

      @exception Exception::InvalidValue is thrown for unregistered names
    */
    void setDescription(const String& name, const String& description);

    /**
      @brief Sets the unit (String), corresponding to an index

      @exception Exception::InvalidValue is thrown for unregistered indices
    */
    void setUnit(UInt index, const String& unit);

    /**
      @brief Sets the unit (String), corresponding to a name

      @exception Exception::InvalidValue is thrown for unregistered names
    */
    void setUnit(const String& name, const String& unit);

    /**
      Returns the integer index corresponding to a string. If the string is not
      registered, returns UInt(-1) (= UINT_MAX).
    */
    UInt getIndex(const String& name) const;

    /**
      @brief Returns the corresponding name to an index

      @exception Exception::InvalidValue is thrown for unregistered indices
    */
    String getName(UInt index) const;

    /**
      @brief returns the description of an index

      @exception Exception::InvalidValue is thrown for unregistered indices
    */
    String getDescription(UInt index) const;
    /**
      @brief returns the description of a name

      @exception Exception::InvalidValue is thrown for unregistered names
    */
    String getDescription(const String& name) const;

    /**
      @brief returns the unit of an index

      @exception Exception::InvalidValue is thrown for unregistered indices
    */
    String getUnit(UInt index) const;
    /**
      @brief returns the unit of a name

      @exception Exception::InvalidValue is thrown for unregistered names
    */
    String getUnit(const String& name) const;

private:
    /// internal counter, that stores the next index to assign
    UInt next_index_;
    using MapString2IndexType = std::unordered_map<std::string, UInt>;
    using MapIndex2StringType = std::unordered_map<UInt, std::string>;
    
    /// map from name to index
    MapString2IndexType name_to_index_;
    /// map from index to name
    MapIndex2StringType index_to_name_;
    /// map from index to description
    MapIndex2StringType index_to_description_;
    /// map from index to unit
    MapIndex2StringType index_to_unit_;
  };

} // namespace OpenMS

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( pop )
#endif

