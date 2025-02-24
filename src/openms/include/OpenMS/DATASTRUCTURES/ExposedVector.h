// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <cstddef>   // for size_t
#include <vector>

namespace OpenMS
{

/// Macro to expose common dependent types, such as @p iterator in the derived class
#define EXPOSED_VECTOR_INTERFACE(InnerElement) \
    using ExpVec = ExposedVector< InnerElement >;                     \
    using ExpVec::ExposedVector;                                               \
    using value_type = typename ExpVec::value_type;                            \
    using iterator = typename ExpVec::iterator;                                \
    using const_iterator = typename ExpVec::const_iterator;                    \
    using reverse_iterator = typename ExpVec::reverse_iterator;                \
    using const_reverse_iterator = typename ExpVec::const_reverse_iterator;    \
    using size_type = typename ExpVec::size_type;                              \
    using pointer = typename ExpVec::pointer;                                  \
    using reference = typename ExpVec::reference;                              \
    using const_reference = typename ExpVec::const_reference;                  \
    using difference_type = typename ExpVec::difference_type;                  \
 

  /**
    @brief Makes a vector<VectorElement> available in the derived class and exposed commonly
           used vector member functions at class level.

    This saves writing repetitive code which forwards commonly used functions of a data member, e.g. 'data_.begin()' 
    as a member function of the class. Also it makes private inheritance from vector<VectorElement> obsolete.
    The latter is problematic for many reasons (read up on 'prefer composition over inheritance'). In our case, even linking
    can be problematic with private inheritance once you require RTTI (which some tools do, e.g. softwipe).

    To fully utilize this class (i.e. access the 'iterator' type), insert
\code
    EXPOSED_VECTOR_INTERFACE(VectorElement)
\endcode

   in your derived class, where @p VectorElement is identical to the template argument of ExposedVector, e.g. 'Feature' for FeatureMap.

   @ingroup Datastructures
  */

  template<class VectorElement>
  class ExposedVector
  {
  public: 
    using VecMember = std::vector<VectorElement>;

    // types
    using value_type = typename VecMember::value_type;
    using iterator = typename VecMember::iterator;
    using const_iterator = typename VecMember::const_iterator;
    using reverse_iterator = typename VecMember::reverse_iterator;
    using const_reverse_iterator = typename VecMember::const_reverse_iterator;
    using size_type = typename VecMember::size_type;
    using pointer = typename VecMember::pointer;
    using reference = typename VecMember::reference;
    using const_reference = typename VecMember::const_reference;
    using difference_type = typename VecMember::difference_type;

  protected:
    VecMember data_; ///< the container which holds all the data

  public:
    ExposedVector() = default;
    explicit ExposedVector(const size_t n)
      : data_(n)
    {
    }
    ExposedVector(const size_t n, const VectorElement& val) : data_(n, val)
    {
    }
    template <class Iter> ExposedVector(Iter begin, Iter end) 
      : data_(begin, end)
    {
    }
    /// Copy C'tor
    ExposedVector(const ExposedVector& rhs) = default;
    /// Move C'tor
    ExposedVector(ExposedVector&& rhs) noexcept = default;

    /// Assignment
    ExposedVector& operator=(const ExposedVector& rhs) = default;
    /// Move Assignment
    ExposedVector& operator=(ExposedVector&& rhs) noexcept = default;

    iterator begin() noexcept
    {
      return data_.begin();
    }
    iterator end() noexcept
    {
      return data_.end();
    }
    const_iterator begin() const noexcept
    {
      return data_.begin();
    }
    const_iterator end() const noexcept
    {
      return data_.end();
    }
    const_iterator cbegin() const noexcept
    {
      return data_.cbegin();
    }
    const_iterator cend() const noexcept
    {
      return data_.cend();
    }
    size_t size() const noexcept
    {
      return data_.size();
    }
    void resize(const size_t new_size)
    {
      data_.resize(new_size);
    }
    void reserve(const size_t new_size)
    {
      data_.reserve(new_size);
    }
    bool empty() const noexcept
    {
      return data_.empty();
    }
    VectorElement& operator[](size_t i) noexcept
    {
      return data_[i];
    }
    const VectorElement& operator[](size_t i) const noexcept
    {
      return data_[i];
    }
    VectorElement& at(size_t i)
    {
      return data_.at(i);
    }
    const VectorElement& at(size_t i) const
    {
      return data_.at(i);
    }
    VectorElement& back() noexcept
    {
      return data_.back();
    }
    const VectorElement& back() const noexcept
    {
      return data_.back();
    }
    void push_back(const VectorElement& f)
    {
      data_.push_back(f);
    }
    void push_back(VectorElement&& f)
    {
      data_.push_back(std::move(f));
    }
    template<typename... Args>
    decltype(auto) emplace_back(Args&&... args)
    {
      return data_.emplace_back(std::forward<Args>(args)...);
    }
    void pop_back() noexcept
    {
      data_.pop_back();
    }
    iterator erase(const_iterator where) noexcept
    {
      return data_.erase(where);
    }
    iterator erase(const_iterator from, const_iterator to) noexcept
    {
      return data_.erase(from, to);
    }
    template<typename T>
    iterator insert(const_iterator where, T from, T to)
    {
      return data_.insert(where, from, to);
    }

    /// read-only access to the underlying data
    const VecMember& getData() const
    {
      return data_;
    }
    /// read access to the underlying data
    VecMember& getData()
    {
      return data_;
    }
  };

} // namespace OpenMS