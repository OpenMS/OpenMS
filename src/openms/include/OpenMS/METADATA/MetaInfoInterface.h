// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

    /// Swap contents
    void swap(MetaInfoInterface& rhs);

    /// Equality operator
    bool operator==(const MetaInfoInterface& rhs) const;
    /// Equality operator
    bool operator!=(const MetaInfoInterface& rhs) const;

    /// Returns the value corresponding to a string, or DataValue::EMPTY if not found
    const DataValue& getMetaValue(const String& name) const;

    /// Returns the value corresponding to a string, or a default value (e.g.: DataValue::EMPTY) if not found    
    DataValue getMetaValue(const String& name, const DataValue& default_value) const; // Note: return needs to be by value to prevent life-time issues at caller site (e.g. if he passes a temporary to default-value)

    /// Returns the value corresponding to the index, or DataValue::EMPTY if not found
    const DataValue& getMetaValue(UInt index) const;

    /// Returns the value corresponding to the index, or a default value (e.g.: DataValue::EMPTY) if not found    
    DataValue getMetaValue(UInt index, const DataValue& default_value) const; // Note: return needs to be by value to prevent life-time issues at caller site

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

    /// function to copy all meta values from one object to this one
    void addMetaValues(const MetaInfoInterface& from);

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
