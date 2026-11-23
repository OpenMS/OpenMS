// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfo.h>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  MetaInfoRegistry MetaInfo::registry_ = MetaInfoRegistry();

  MetaInfo::iterator MetaInfo::find_(UInt index)
  {
    auto it = std::lower_bound(begin(), end(), index, [](const auto& entry, UInt key) { return entry.first < key; });
    return it != end() && it->first == index ? it : end();
  }

  MetaInfo::const_iterator MetaInfo::find_(UInt index) const
  {
    auto it = std::lower_bound(begin(), end(), index, [](const auto& entry, UInt key) { return entry.first < key; });
    return it != end() && it->first == index ? it : end();
  }

  MetaInfo::~MetaInfo() = default;

  bool MetaInfo::operator==(const MetaInfo& rhs) const
  {
    return index_to_value_ == rhs.index_to_value_;
  }

  bool MetaInfo::operator!=(const MetaInfo& rhs) const
  {
    return !(operator==(rhs));
  }

  MetaInfo& MetaInfo::operator+=(const MetaInfo& rhs)
  {
    if (rhs.index_to_value_.empty()) return *this;
    if (index_to_value_.empty())
    {
      index_to_value_ = rhs.index_to_value_;
      return *this;
    }

    // Two-way merge into a sorted vector
    using pair_type = MapType::value_type;
    std::vector<pair_type> merged;
    merged.reserve(index_to_value_.size() + rhs.index_to_value_.size());

    auto it_this = index_to_value_.begin();
    auto end_this = index_to_value_.end();
    auto it_rhs = rhs.index_to_value_.begin();
    auto end_rhs = rhs.index_to_value_.end();

    // Merge with deduplication: rhs values overwrite
    while (it_this != end_this && it_rhs != end_rhs)
    {
      if (it_this->first < it_rhs->first)
      {
        merged.push_back(*it_this++);
      }
      else if (it_rhs->first < it_this->first)
      {
        merged.push_back(*it_rhs++);
      }
      else
      {
        // Equal keys: rhs value overwrites
        merged.push_back(*it_rhs++);
        ++it_this;
      }
    }
    // Append remaining elements efficiently
    merged.insert(merged.end(), it_this, end_this);
    merged.insert(merged.end(), it_rhs, end_rhs);

    // Keep the merged storage without another element-wise copy
    index_to_value_ = std::move(merged);
    return *this;
  }

  const DataValue& MetaInfo::getValue(const std::string& name, const DataValue& default_value) const
  {
    MapType::const_iterator it = find_(registry_.getIndex(name));
    if (it != index_to_value_.end())
    {
      return it->second;
    }
    return default_value;
  }

  const DataValue& MetaInfo::getValue(UInt index, const DataValue& default_value) const
  {
    MapType::const_iterator it = find_(index);
    if (it != index_to_value_.end())
    {
      return it->second;
    }
    return default_value;
  }

  void MetaInfo::setValue(const std::string& name, const DataValue& value)
  {
    UInt index = registry_.registerName(name); // no-op if name is already registered
    setValue(index, value);
  }

  void MetaInfo::setValue(UInt index, const DataValue& value)
  {
    // @TODO: check if that index is registered in MetaInfoRegistry?
    auto it = find_(index);
    if (it != index_to_value_.end())
    {
      it->second = value;
    }
    else
    {
      // Note; we need to create a copy of data value here and can't use the const &
      // The underlying sorted vector invalidates references to it if inserting
      // an element leads to relocation (e.g, in constructs like: m.insert(1, m[2]));)
      DataValue tmp = value;
      auto pos
        = std::lower_bound(index_to_value_.begin(), index_to_value_.end(), index, [](const auto& entry, UInt key) { return entry.first < key; });
      index_to_value_.insert(pos, std::make_pair(index, std::move(tmp)));
    }
  }

  MetaInfoRegistry& MetaInfo::registry()
  {
    return registry_;
  }

  bool MetaInfo::exists(const std::string& name) const
  {
    UInt index = registry_.getIndex(name);
    if (index != UInt(-1)) { return (find_(index) != index_to_value_.end()); }
    return false;
  }

  bool MetaInfo::exists(UInt index) const
  { return (find_(index) != index_to_value_.end()); }

  void MetaInfo::removeValue(const std::string& name)
  {
    MapType::iterator it = find_(registry_.getIndex(name));
    if (it != index_to_value_.end())
    {
      index_to_value_.erase(it);
    }
  }

  void MetaInfo::removeValue(UInt index)
  {
    MapType::iterator it = find_(index);
    if (it != index_to_value_.end())
    {
      index_to_value_.erase(it);
    }
  }

  void MetaInfo::getKeys(vector<std::string>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i = 0;
    for (MapType::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
    {
      keys[i++] = registry_.getName(it->first);
    }
  }

  void MetaInfo::getKeys(vector<UInt>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i = 0;
    for (MapType::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
    {
      keys[i++] = it->first;
    }
  }

  bool MetaInfo::empty() const
  {
    return index_to_value_.empty();
  }

  void MetaInfo::clear()
  {
    index_to_value_.clear();
  }

} //namespace

