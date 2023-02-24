// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <vector>

namespace OpenMS
{
  /**
    @brief Makes a vector<VectorElement> available in the derived class and exposed commonly
           used vector member functions at class level.

    This saves writing repetitive code which exposes 'data_.begin()' as a member function 
    and makes private inheritance from vector<VectorElement> obsolete. The latter is problematic
    for many reasons (read up on 'prefer composition over inheritance'). In our case, even linking
    can be problematic with private inheritance once you require RTTI (which some tools do, e.g. softwipe).


    @ingroup Datastructures
  */

  template<class VectorElement>
  class ExposedVector
  {
    template<class T> friend class UniqueIdIndexer; ///< accesses data_ but is not derived from ExposedVector

    using Base = std::vector<VectorElement>;

  protected:
    Base data_;  ///< the container which holds all the data

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

    typename Base::iterator begin() noexcept
    {
      return data_.begin();
    }
    typename Base::iterator end() noexcept
    {
      return data_.end();
    }
    typename Base::const_iterator begin() const noexcept
    {
      return data_.begin();
    }
    typename Base::const_iterator end() const noexcept
    {
      return data_.end();
    }
    typename Base::const_iterator cbegin() const noexcept
    {
      return data_.cbegin();
    }
    typename Base::const_iterator cend() const noexcept
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
    typename Base::iterator erase(typename Base::const_iterator where) noexcept
    {
      return data_.erase(where);
    }
    typename Base::iterator erase(typename Base::const_iterator from, typename Base::const_iterator to) noexcept
    {
      return data_.erase(from, to);
    }
    template<typename T>
    typename Base::iterator insert(typename Base::const_iterator where, T from, T to)
    {
      return data_.insert(where, from, to);
    }

    // types
    using value_type = typename Base::value_type;
    using iterator = typename Base::iterator;
    using const_iterator = typename Base::const_iterator;
    using size_type = typename Base::size_type;
    using pointer = typename Base::pointer;                 // ConstRefVector
    using reference = typename Base::reference;             // ConstRefVector
    using const_reference = typename Base::const_reference; // ConstRefVector
    using difference_type = typename Base::difference_type; // ConstRefVector
  };

} // namespace OpenMS