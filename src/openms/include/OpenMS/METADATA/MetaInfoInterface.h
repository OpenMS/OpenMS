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

namespace OpenMS
{
  class String;
  class MetaInfo;

  /**
    @brief Interface for classes that can store arbitrary meta information
    (Type-Name-Value tuples).

    MetaInfoInterface is a base class for all classes that use one MetaInfo
    object as member.  If you want to add meta information to a class, let it
    publicly inherit the MetaInfoInterface.  Meta information is an array of
    Type-Name-Value tuples.

    @ingroup Metadata
  */

  class OPENMS_DLLAPI MetaInfoInterface
  {
public:

    /// Constructor
    MetaInfoInterface();
    /// Copy constructor
    MetaInfoInterface(const MetaInfoInterface& rhs);
    /// Move constructor
    MetaInfoInterface(MetaInfoInterface&&) noexcept;
    /// Destructor
    ~MetaInfoInterface();

    /// Assignment operator
    MetaInfoInterface& operator=(const MetaInfoInterface& rhs);
    /// Move assignment operator
    MetaInfoInterface& operator=(MetaInfoInterface&&) noexcept;

    /// Equality operator
    bool operator==(const MetaInfoInterface& rhs) const;
    /// Equality operator
    bool operator!=(const MetaInfoInterface& rhs) const;

    /// Returns the value corresponding to a string, or a default value (default: DataValue::EMPTY) if not found
    const DataValue& getMetaValue(const String& name, const DataValue& default_value = DataValue::EMPTY) const;
    /// Returns the value corresponding to an index, or a default value (default: DataValue::EMPTY) if not found
    const DataValue& getMetaValue(UInt index, const DataValue& default_value = DataValue::EMPTY) const;

    /// Returns whether an entry with the given name exists
    bool metaValueExists(const String& name) const;
    /// Returns whether an entry with the given index exists
    bool metaValueExists(UInt index) const;

    /// Sets the DataValue corresponding to a name
    void setMetaValue(const String& name, const DataValue& value);
    /// Sets the DataValue corresponding to an index
    void setMetaValue(UInt index, const DataValue& value);

    /// Removes the DataValue corresponding to @p name if it exists
    void removeMetaValue(const String& name);
    /// Removes the DataValue corresponding to @p index if it exists
    void removeMetaValue(UInt index);

    /// Returns a reference to the MetaInfoRegistry
    static MetaInfoRegistry& metaRegistry();

    /// Fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<String>& keys) const;

    /// Fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<UInt>& keys) const;

    /// Returns if the MetaInfo is empty
    bool isMetaEmpty() const;

    /// Removes all meta values
    void clearMetaInfo();

protected:

    /// Creates the MetaInfo object if it does not exist
    inline void createIfNotExists_();

    /// Pointer to the MetaInfo object (0 by default)
    MetaInfo* meta_;
  };

} // namespace OpenMS

