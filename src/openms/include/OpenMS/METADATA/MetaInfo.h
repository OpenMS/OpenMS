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

#include <boost/container/flat_map.hpp>
#include <functional>

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
    /// Internal map type (UInt key to DataValue)
    using MapType = boost::container::flat_map<UInt, DataValue>;
    /// Mutable iterator type
    using iterator = MapType::iterator;
    /// Const iterator type
    using const_iterator = MapType::const_iterator;

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

    /**
     * @brief Merge another MetaInfo into this one.
     *
     * All entries from @p rhs are added to this MetaInfo.
     * If an entry with the same index already exists, it will be overwritten
     * with the value from @p rhs.
     *
     * Uses an O(n+m) two-way merge algorithm since the underlying flat_map is sorted.
     *
     * @param rhs The MetaInfo to merge from.
     * @return Reference to this object.
     */
    MetaInfo& operator+=(const MetaInfo& rhs);

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

    /// @name Iterator access
    /// @brief Provides iterator access to the underlying index-to-value map.
    /// Iterators dereference to std::pair<UInt, DataValue> where the first element
    /// is the registry index and the second is the associated value.
    /// The iteration order is sorted by index (ascending).
    ///@{

    /**
     * @brief Returns a const iterator to the beginning of the meta info entries.
     * @return const_iterator pointing to the first index-value pair, or end() if empty.
     */
    const_iterator begin() const { return index_to_value_.begin(); }

    /**
     * @brief Returns a const iterator to the end of the meta info entries.
     * @return const_iterator pointing past the last index-value pair.
     */
    const_iterator end() const { return index_to_value_.end(); }

    /**
     * @brief Returns a const iterator to the beginning of the meta info entries.
     * @return const_iterator pointing to the first index-value pair, or cend() if empty.
     */
    const_iterator cbegin() const { return index_to_value_.cbegin(); }

    /**
     * @brief Returns a const iterator to the end of the meta info entries.
     * @return const_iterator pointing past the last index-value pair.
     */
    const_iterator cend() const { return index_to_value_.cend(); }

    /**
     * @brief Returns a mutable iterator to the beginning of the meta info entries.
     * @return iterator pointing to the first index-value pair, or end() if empty.
     * @note Modifying the key (first element) may invalidate the sorted order invariant.
     */
    iterator begin() { return index_to_value_.begin(); }

    /**
     * @brief Returns a mutable iterator to the end of the meta info entries.
     * @return iterator pointing past the last index-value pair.
     */
    iterator end() { return index_to_value_.end(); }
    ///@}

    /**
     * @brief Returns the number of meta value entries.
     * @return The count of index-value pairs stored.
     */
    Size size() const { return index_to_value_.size(); }

private:

    /// Static MetaInfoRegistry
    static MetaInfoRegistry registry_;

    /// The actual mapping of indexes to values
    MapType index_to_value_;

    // Grant access to hash implementation
    friend struct std::hash<MetaInfo>;
  };

} // namespace OpenMS

// Hash function specialization for MetaInfo
namespace std
{
  template<>
  struct hash<OpenMS::MetaInfo>
  {
    std::size_t operator()(const OpenMS::MetaInfo& mi) const noexcept
    {
      std::size_t seed = 0;
      // Hash each key-value pair in sorted order (flat_map is already sorted by key)
      for (const auto& [key, value] : mi.index_to_value_)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_int(key));
        OpenMS::hash_combine(seed, std::hash<OpenMS::DataValue>{}(value));
      }
      return seed;
    }
  };
} // namespace std

