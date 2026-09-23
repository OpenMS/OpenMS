// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

namespace OpenMS::IdentificationDataInternal
{
/**
  @brief An indexed collection of identification records with stable references.

  The implementation retains ordered uniqueness, constant-time iterator movement,
  and stable references on insertion and successful modification. The underlying
  indexing library is private to OpenMS. A modification that violates uniqueness
  erases the modified record and returns false, invalidating references to it.
  Copying copies records; moving and swapping preserve references to records.

  Only the explicitly instantiated identification record types are supported.
*/
template<typename Value, typename Key, typename Prefix = Key>
class IDDataContainer
{
  struct Impl;

public:
  using value_type = Value;
  using key_type = Key;
  using size_type = std::size_t;

  class iterator
  {
    friend class IDDataContainer;

  protected:
    const Impl* owner_ = nullptr;
    const Value* value_ = nullptr;
    iterator(const Impl* owner, const Value* value): owner_(owner), value_(value)
    {
    }

  public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = Value;
    using difference_type = std::ptrdiff_t;
    using pointer = const Value*;
    using reference = const Value&;
    iterator() = default;
    reference operator*() const
    { return *value_; }
    pointer operator->() const
    { return value_; }
    iterator& operator++()
    {
      value_ = advance_(owner_, value_, true, false);
      return *this;
    }
    iterator operator++(int)
    {
      auto old = *this;
      ++*this;
      return old;
    }
    iterator& operator--()
    {
      value_ = advance_(owner_, value_, false, false);
      return *this;
    }
    iterator operator--(int)
    {
      auto old = *this;
      --*this;
      return old;
    }
    bool operator==(const iterator& other) const
    { return owner_ == other.owner_ && value_ == other.value_; }
  };
  /// Iterator for the ordered processing-step view, with no per-iterator allocation.
  class ordered_iterator : public iterator
  {
    friend class IDDataContainer;
    ordered_iterator(const Impl* owner, const Value* value): iterator(owner, value)
    {
    }

  public:
    ordered_iterator() = default;
    ordered_iterator& operator++()
    {
      this->value_ = advance_(this->owner_, this->value_, true, true);
      return *this;
    }
    ordered_iterator operator++(int)
    {
      auto old = *this;
      ++*this;
      return old;
    }
    ordered_iterator& operator--()
    {
      this->value_ = advance_(this->owner_, this->value_, false, true);
      return *this;
    }
    ordered_iterator operator--(int)
    {
      auto old = *this;
      --*this;
      return old;
    }
  };
  using const_iterator = iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = reverse_iterator;

  IDDataContainer();
  IDDataContainer(const IDDataContainer& other);
  IDDataContainer(IDDataContainer&& other) noexcept;
  IDDataContainer& operator=(const IDDataContainer& other);
  IDDataContainer& operator=(IDDataContainer&& other) noexcept;
  ~IDDataContainer();
  IDDataContainer(std::initializer_list<Value> values): IDDataContainer()
  {
    for (const auto& value : values)
      insert(value);
  }

  iterator begin() const;
  iterator end() const
  { return iterator(impl_.get(), nullptr); }
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
  { return size() == 0; }
  size_type size() const;
  void clear();
  void swap(IDDataContainer& other) noexcept
  { impl_.swap(other.impl_); }
  std::pair<iterator, bool> insert(const Value& value);
  std::pair<iterator, bool> push_back(const Value& value)
  { return insert(value); }
  template<typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args)
  { return insert(Value(std::forward<Args>(args)...)); }
  template<typename... Args>
  std::pair<iterator, bool> emplace_back(Args&&... args)
  { return emplace(std::forward<Args>(args)...); }
  iterator erase(iterator position);
  /// Changes the record at @p position in place by calling @p modifier(Value&); see the class
  /// description for a change that breaks uniqueness.
  template<typename Modifier>
  bool modify(iterator position, Modifier&& modifier)
  {
    // Passed to the out-of-line modify_() as a plain function pointer and a pointer to the
    // caller's callable: unlike std::function, this never copies or allocates the callable.
    using Callable = std::remove_reference_t<Modifier>;
    auto call = [](void* context, Value& value) { (*static_cast<Callable*>(context))(value); };
    return modify_(position, call, const_cast<void*>(static_cast<const void*>(std::addressof(modifier))));
  }
  iterator find(const Key& key) const;
  std::pair<iterator, iterator> equal_range(const Prefix& key) const;

  /// Ordered view of a sequenced collection (applied processing steps).
  class ConstOrderedView
  {
  protected:
    const IDDataContainer* owner_;

  public:
    explicit ConstOrderedView(const IDDataContainer* owner): owner_(owner)
    {
    }
    ordered_iterator begin() const
    { return owner_->orderedBegin_(); }
    ordered_iterator end() const
    { return ordered_iterator(owner_->impl_.get(), nullptr); }
    ordered_iterator find(const Key& key) const
    {
      auto position = owner_->find(key);
      return ordered_iterator(owner_->impl_.get(), position == owner_->end() ? nullptr : std::addressof(*position));
    }
    size_type size() const
    { return owner_->size(); }
    bool empty() const
    { return owner_->empty(); }
  };
  class OrderedView : public ConstOrderedView
  {
    IDDataContainer* mutable_owner_;

  public:
    explicit OrderedView(IDDataContainer* owner): ConstOrderedView(owner), mutable_owner_(owner)
    {
    }
    template<typename Modifier>
    bool modify(iterator position, Modifier&& modifier)
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
  friend bool operator==(const IDDataContainer& lhs, const IDDataContainer& rhs)
    requires std::equality_comparable<Value>
  { return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin()); }

private:
  std::unique_ptr<Impl> impl_;
  ordered_iterator orderedBegin_() const;
  bool modify_(iterator position, void (*call)(void*, Value&), void* context);
  static const Value* advance_(const Impl* owner, const Value* value, bool forward, bool ordered);
};
} // namespace OpenMS::IdentificationDataInternal
