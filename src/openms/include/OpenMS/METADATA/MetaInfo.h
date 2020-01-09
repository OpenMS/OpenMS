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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <boost/container/flat_map.hpp>

namespace OpenMS
{
  class String;

  /**
      @brief A Type-Name-Value tuple class.

      MetaInfo maps an index (an integer corresponding to a string) to
      DataValue objects.  The mapping of strings to the index is performed by
      the MetaInfoRegistry, which can be accessed by the method registry().

      There are two versions of nearly all members. One which operates with a
      string name and another one which operates on an index. The index version
      is always faster, as it does not need to look up the index corresponding
      to the string in the MetaInfoRegistry.

      If you wish to add a MetaInfo member to a class, consider deriving that
      class from MetaInfoInterface, instead of simply adding MetaInfo as
      member. MetaInfoInterface implements a full interface to a MetaInfo
      member and is more memory efficient if no meta info gets added.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI MetaInfo
  {
public:
    /// Constructor
    MetaInfo() = default;

    /// Copy constructor
    MetaInfo(const MetaInfo&) = default;

    /// Move constructor
    MetaInfo(MetaInfo&&) = default;

    /// Destructor
    ~MetaInfo();

    /// Assignment operator
    MetaInfo& operator=(const MetaInfo&) = default;
    /// Move assignment operator
    MetaInfo& operator=(MetaInfo&&) & = default;

    /// Equality operator
    bool operator==(const MetaInfo& rhs) const;
    /// Equality operator
    bool operator!=(const MetaInfo& rhs) const;

    /// Returns the value corresponding to a string, or a default value (default: DataValue::EMPTY) if not found
    const DataValue& getValue(const String& name, const DataValue& default_value = DataValue::EMPTY) const;
    /// Returns the value corresponding to an index, or a default value (default: DataValue::EMPTY) if not found
    const DataValue& getValue(UInt index, const DataValue& default_value = DataValue::EMPTY) const;

    /// Returns whether an entry with the given name exists
    bool exists(const String& name) const;
    /// Returns whether an entry with the given index exists
    bool exists(UInt index) const;

    /// Sets the DataValue corresponding to a name
    void setValue(const String& name, const DataValue& value);
    /// Sets the DataValue corresponding to an index
    void setValue(UInt index, const DataValue& value);

    /// Removes the DataValue corresponding to @p name if it exists
    void removeValue(const String& name);
    /// Removes the DataValue corresponding to @p index if it exists
    void removeValue(UInt index);

    /// Returns a reference to the MetaInfoRegistry
    static MetaInfoRegistry& registry();

    /// Fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<String>& keys) const;

    /// Fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<UInt>& keys) const;

    /// Returns if the MetaInfo is empty
    bool empty() const;

    /// Removes all meta values
    void clear();

private:
    using MapType = boost::container::flat_map<UInt, DataValue>;

    /// Static MetaInfoRegistry
    static MetaInfoRegistry registry_;

    /// The actual mapping of indexes to values
    MapType index_to_value_;
  };

} // namespace OpenMS

