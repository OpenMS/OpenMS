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
#include <OpenMS/CONCEPT/HashUtils.h>
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
    MetaInfoInterface() = default;
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

    /// Copy all meta values from @p from to this object (both named String and indexed UInt keys).
    /// Existing entries with the same key are overwritten; others are preserved.
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

    /// Pointer to the MetaInfo object
    MetaInfo* meta_ = nullptr;
  };

} // namespace OpenMS

// Hash function specialization for MetaInfoInterface
namespace std
{
  /**
   * @brief Hash function for MetaInfoInterface.
   *
   * Hashes based on all meta info key-value pairs using order-independent
   * accumulation. This handles both String keys (names) and UInt keys (indices).
   *
   * @note Uses additive accumulation (commutative) instead of sorting for O(n)
   *       complexity. Each key-value pair produces a hash that is added to the
   *       total, ensuring equal MetaInfoInterface objects produce equal hashes
   *       regardless of internal storage order.
   */
  template<>
  struct hash<OpenMS::MetaInfoInterface>
  {
    std::size_t operator()(const OpenMS::MetaInfoInterface& meta) const noexcept
    {
      // Use additive accumulation for order-independent hashing (O(n) vs O(n log n))
      std::size_t hash = 0;

      // Hash String keys (names)
      std::vector<OpenMS::String> str_keys;
      meta.getKeys(str_keys);
      for (const auto& key : str_keys)
      {
        std::size_t pair_hash = OpenMS::fnv1a_hash_string(key);
        OpenMS::hash_combine(pair_hash, std::hash<OpenMS::DataValue>{}(meta.getMetaValue(key)));
        hash += pair_hash;  // Order-independent accumulation
      }

      // Hash UInt keys (indices) - critical for correctness
      std::vector<OpenMS::UInt> uint_keys;
      meta.getKeys(uint_keys);
      for (const auto& key : uint_keys)
      {
        std::size_t pair_hash = OpenMS::hash_int(static_cast<int64_t>(key));
        OpenMS::hash_combine(pair_hash, std::hash<OpenMS::DataValue>{}(meta.getMetaValue(key)));
        hash += pair_hash;  // Order-independent accumulation
      }

      return hash;
    }
  };
} // namespace std
