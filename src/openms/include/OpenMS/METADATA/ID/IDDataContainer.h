// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace OpenMS::IdentificationDataInternal
{
namespace Detail
{
  template<typename T>
  inline constexpr bool is_tuple = false;
  template<typename... T>
  inline constexpr bool is_tuple<std::tuple<T...>> = true;

  /// Lexicographic "less" over the leading elements two tuples have in common, so a key
  /// prefix (e.g. the input file of an observation) compares against a full key.
  template<std::size_t I = 0, typename A, typename B>
  bool tupleLess(const A& a, const B& b)
  {
    constexpr std::size_t n = std::min(std::tuple_size_v<A>, std::tuple_size_v<B>);
    if constexpr (I == n) return false;
    else
    {
      if (std::get<I>(a) < std::get<I>(b)) return true;
      if (std::get<I>(b) < std::get<I>(a)) return false;
      return tupleLess<I + 1>(a, b);
    }
  }
} // namespace Detail

/**
  @brief A collection of identification records, unique and ordered by the key formed from
  the members @p KeyMembers.

  A std::set underneath: references and iterators to records stay valid on insertion and on a
  modification that keeps the record's key, and they are 8-byte iterators that step inline.
  Lookup takes the key (a std::tuple for a composite key) or a @p Prefix of it. Copying copies
  the records; moving and swapping keep references to them.

  modify() changes a record in place. A modification that changes the key moves the record to
  its new position (references to it stay valid); one that duplicates another record's key
  erases the record instead and returns false, as Boost.MultiIndex did.
*/
template<typename Value, typename Key, typename Prefix = Key, auto... KeyMembers>
class IDDataContainer
{
  static_assert(sizeof...(KeyMembers) > 0, "the key members of the record type are required");

  // A record compares by its key members, a key or prefix by itself, both as tuples.
  template<typename T>
  static decltype(auto) keyView(const T& x)
  {
    if constexpr (std::is_same_v<T, Value>) return std::tie(x.*KeyMembers...);
    else if constexpr (Detail::is_tuple<T>) return (x);
    else return std::tie(x);
  }
  struct Compare
  {
    using is_transparent = void;
    template<typename A, typename B>
    bool operator()(const A& a, const B& b) const
    { return Detail::tupleLess(keyView(a), keyView(b)); }
  };
  using Set = std::set<Value, Compare>;
  Set records_;

public:
  using value_type = Value;
  using key_type = Key;
  using size_type = std::size_t;
  using iterator = typename Set::const_iterator;
  using const_iterator = iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = reverse_iterator;

  IDDataContainer() = default;
  IDDataContainer(std::initializer_list<Value> values)
  {
    for (const auto& value : values)
      insert(value);
  }

  iterator begin() const
  { return records_.begin(); }
  iterator end() const
  { return records_.end(); }
  iterator cbegin() const
  { return begin(); }
  iterator cend() const
  { return end(); }
  reverse_iterator rbegin() const
  { return reverse_iterator(end()); }
  reverse_iterator rend() const
  { return reverse_iterator(begin()); }
  const Value& front() const
  { return *begin(); }
  const Value& back() const
  { return *--end(); }
  bool empty() const
  { return records_.empty(); }
  size_type size() const
  { return records_.size(); }
  void clear()
  { records_.clear(); }
  void swap(IDDataContainer& other) noexcept
  { records_.swap(other.records_); }
  std::pair<iterator, bool> insert(const Value& value)
  { return records_.insert(value); }
  std::pair<iterator, bool> push_back(const Value& value)
  { return insert(value); }
  template<typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args)
  { return insert(Value(std::forward<Args>(args)...)); }
  template<typename... Args>
  std::pair<iterator, bool> emplace_back(Args&&... args)
  { return emplace(std::forward<Args>(args)...); }
  iterator erase(iterator position)
  { return records_.erase(position); }
  /**
    @brief Calls @p modifier(Value&) on the record at @p position.

    @p position follows the record: it keeps pointing at it after a key change, and after a
    modification that was rejected as a duplicate it points at the record after the erased one.
  */
  template<typename Modifier>
  bool modify(iterator& position, Modifier&& modifier)
  {
    // The set exposes its records as const to protect their order; the key is checked below.
    modifier(const_cast<Value&>(*position));
    if (stillOrdered_(position)) return true;
    auto next = std::next(position);
    auto node = records_.extract(position);
    auto result = records_.insert(std::move(node));
    if (result.inserted)
    {
      position = result.position;
      return true;
    }
    position = next;
    return false;
  }
  iterator find(const Key& key) const
  { return records_.find(key); }
  std::pair<iterator, iterator> equal_range(const Prefix& key) const
  { return records_.equal_range(key); }

  /// Value equality is available only for equality-comparable record types.
  friend bool operator==(const IDDataContainer& lhs, const IDDataContainer& rhs)
    requires std::equality_comparable<Value>
  { return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin()); }

private:
  /// True if the record at @p position still sorts strictly between its neighbours.
  bool stillOrdered_(iterator position) const
  {
    const auto& less = records_.key_comp();
    if (position != records_.begin() && ! less(*std::prev(position), *position)) return false;
    auto next = std::next(position);
    return next == records_.end() || less(*position, *next);
  }
};

/**
  @brief Records kept in the order they were added, unique by the member @p KeyMember.

  Used for the processing steps applied to a record: a handful per record, looked up by step.
  get<1>() gives the lookup view by key that Boost.MultiIndex's second index used to provide.
*/
template<typename Value, typename Key, auto KeyMember>
class IDSequencedContainer
{
  std::vector<Value> records_;

public:
  using value_type = Value;
  using key_type = Key;
  using size_type = std::size_t;
  using iterator = typename std::vector<Value>::const_iterator;
  using const_iterator = iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = reverse_iterator;

  IDSequencedContainer() = default;
  IDSequencedContainer(std::initializer_list<Value> values)
  {
    for (const auto& value : values)
      push_back(value);
  }

  iterator begin() const
  { return records_.begin(); }
  iterator end() const
  { return records_.end(); }
  iterator cbegin() const
  { return begin(); }
  iterator cend() const
  { return end(); }
  reverse_iterator rbegin() const
  { return reverse_iterator(end()); }
  reverse_iterator rend() const
  { return reverse_iterator(begin()); }
  const Value& front() const
  { return records_.front(); }
  const Value& back() const
  { return records_.back(); }
  bool empty() const
  { return records_.empty(); }
  size_type size() const
  { return records_.size(); }
  void clear()
  { records_.clear(); }
  void swap(IDSequencedContainer& other) noexcept
  { records_.swap(other.records_); }
  /// Appends @p value unless a record with its key exists; then that record and false
  std::pair<iterator, bool> push_back(const Value& value)
  {
    auto position = find(value.*KeyMember);
    if (position != end()) return {position, false};
    records_.push_back(value);
    return {std::prev(end()), true};
  }
  std::pair<iterator, bool> insert(const Value& value)
  { return push_back(value); }
  template<typename... Args>
  std::pair<iterator, bool> emplace_back(Args&&... args)
  { return push_back(Value(std::forward<Args>(args)...)); }
  template<typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args)
  { return emplace_back(std::forward<Args>(args)...); }
  iterator erase(iterator position)
  { return records_.erase(position); }
  iterator find(const Key& key) const
  {
    return std::find_if(records_.begin(), records_.end(), [&key](const Value& value) { return value.*KeyMember == key; });
  }
  std::pair<iterator, iterator> equal_range(const Key& key) const
  {
    auto position = find(key);
    return {position, position == end() ? position : std::next(position)};
  }
  /// As IDDataContainer::modify(): a modification that duplicates another record's key erases the record.
  template<typename Modifier>
  bool modify(iterator& position, Modifier&& modifier)
  {
    Value& value = records_[position - begin()];
    modifier(value);
    for (auto it = begin(); it != end(); ++it)
    {
      if (it != position && (*it).*KeyMember == value.*KeyMember)
      {
        position = erase(position);
        return false;
      }
    }
    return true;
  }

  /// Lookup of records by key
  class ConstOrderedView
  {
  protected:
    const IDSequencedContainer* owner_;

  public:
    explicit ConstOrderedView(const IDSequencedContainer* owner): owner_(owner)
    {
    }
    iterator find(const Key& key) const
    { return owner_->find(key); }
    iterator end() const
    { return owner_->end(); }
    size_type size() const
    { return owner_->size(); }
    bool empty() const
    { return owner_->empty(); }
  };
  class OrderedView : public ConstOrderedView
  {
    IDSequencedContainer* mutable_owner_;

  public:
    explicit OrderedView(IDSequencedContainer* owner): ConstOrderedView(owner), mutable_owner_(owner)
    {
    }
    template<typename Modifier>
    bool modify(iterator& position, Modifier&& modifier)
    { return mutable_owner_->modify(position, std::forward<Modifier>(modifier)); }
  };
  template<int Index>
  struct nth_index
  {
    static_assert(Index == 1);
    using type = OrderedView;
  };
  template<int Index>
  OrderedView get()
  {
    static_assert(Index == 1);
    return OrderedView(this);
  }
  template<int Index>
  ConstOrderedView get() const
  {
    static_assert(Index == 1);
    return ConstOrderedView(this);
  }

  /// Value equality is available only for equality-comparable record types.
  friend bool operator==(const IDSequencedContainer& lhs, const IDSequencedContainer& rhs)
    requires std::equality_comparable<Value>
  { return lhs.records_ == rhs.records_; }
};
} // namespace OpenMS::IdentificationDataInternal
