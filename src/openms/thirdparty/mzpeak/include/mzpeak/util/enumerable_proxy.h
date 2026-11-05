/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <boost/range/detail/common.hpp>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <utility>

#include "mzpeak/exception.h"

namespace MzPeak::Util {

template <typename T, typename V = T>
class EnumerableProxy : public std::ranges::view_interface<EnumerableProxy<T>> {
public:
  /// The iterator type.
  class Iterator {
  public:
    using difference_type = std::ptrdiff_t;
    using value_type = V;

    /// Constructor for a valid iterator.
    Iterator(std::size_t n, EnumerableProxy<T, V>* owner)
        : n_(n)
        , owner_(owner)
        , cache_()
    {
    }

    /// Constructor for an invalid iterator.
    Iterator(std::size_t n)
        : n_(n)
        , owner_(nullptr)
        , cache_()
    {
    }

    /// Default construction also invalid;
    Iterator()
        : n_(0)
        , owner_(nullptr)
        , cache_()
    {
    }

    // Copy, move, and assignment constructors.
    //
    // The move overloads take `Iterator&&`, NOT `const Iterator&&`: a defaulted
    // move constructor or move assignment may only take a non-const rvalue
    // reference, so `= default` on the const&& forms is ill-formed.  Clang
    // accepts it, GCC rejects it, so the const versions made this header fail
    // to compile on GCC as soon as anything instantiated the template.
    Iterator(const Iterator&) = default;
    Iterator(Iterator&&) = default;
    Iterator& operator=(const Iterator&) = default;
    Iterator& operator=(Iterator&&) = default;

    /// Prefix increment.
    Iterator& operator++()
    {
      ++n_;
      return *this;
    }

    /// Postfix increment.
    Iterator operator++(int)
    {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    /// Prefix decrement.
    Iterator& operator--()
    {
      --n_;
      return *this;
    }

    /// Postfix decrement.
    Iterator operator--(int)
    {
      auto tmp = *this;
      --*this;
      return tmp;
    }

    /// Equality operator.
    bool operator==(const Iterator& other) const { return n_ == other.n_; }

    /// Dereference operator.
    value_type operator*() const
    {
      if (cache_.has_value() && cache_->first == n_) {
        return cache_->second;
      } else if (owner_ != nullptr) {
        auto* self = const_cast<Iterator*>(this);
        self->cache_ = std::make_pair<>(n_, self->owner_->fetch(n_));
        return cache_->second;
      } else {
        throw InvalidIteratorError("attempt to dereference an invalid iterator");
      }
    }

  private:
    std::size_t n_;
    EnumerableProxy<T, V>* owner_;
    std::optional<std::pair<std::size_t, value_type>> cache_;
  };

  // Ensure the iterator is correctly defined.
  //
  // NOTE: This *could* be a `random_access_iterator` but I'd need to
  // implement six more methods.  If we need it then I'll do it.
  static_assert(std::bidirectional_iterator<Iterator>);

  /// Default constructor.
  EnumerableProxy() = default;

  /// Destructor.
  ~EnumerableProxy() = default;

  /// Meaningful constructor.
  EnumerableProxy(std::size_t count)
      : count_(count)
  {
  }

  /// The number of elements.
  std::size_t size() const { return count_; }

  /// Iterator to the first element.
  Iterator begin() const
  {
    return Iterator(0, const_cast<EnumerableProxy<T, V>*>(this));
  }

  /// Iterator/sentinel representing one element beyond the last.
  Iterator end() const { return Iterator(count_); }

  // Access an element by its index.
  V operator[](std::size_t n) { return fetch(n); }

protected:
  friend class Iterator;

  /// Update the internal count of records.
  void resize(std::size_t n) { count_ = n; }

  /// The function to call to get a specific element.
  virtual V fetch(uint64_t) = 0;

private:
  std::size_t count_ = 0;
};

} // namespace MzPeak::Util
