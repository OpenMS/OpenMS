// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MAP_H
#define OPENMS_DATASTRUCTURES_MAP_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>

#include <map>

namespace OpenMS
{
  /**
    @brief Map class based on the STL map (containing several convenience functions)

    @ingroup Datastructures
  */
  template <class Key, class T>
  class Map :
    private std::map<Key, T>
  {
public:

    /**
      @brief Map illegal key exception

      Thrown when trying to access an element with operator[], which is not contained in the Map
      , i.e. no default ctor is called (as done in std::map)

      @ingroup Exceptions
    */
    class IllegalKey :
      public Exception::BaseException
    {
public:
      IllegalKey(const char* file, int line, const char* function) :
        Exception::BaseException(file, line, function)
      {
      }

    };

    ///@name OpenMS style typedefs
    //@{
    typedef std::map<Key, T> Base;
    typedef typename Base::value_type ValueType;
    typedef Key KeyType;
    typedef typename Base::value_type* PointerType;
    typedef typename Base::iterator Iterator;
    typedef typename Base::const_iterator ConstIterator;
    typedef typename Base::reverse_iterator ReverseIterator;
    typedef typename Base::const_reverse_iterator ConstReverseIterator;
    //@}

    ///@name Export methods from private base std::vector<Acquisition>
    //@{
    using Base::erase;
    using Base::size;
    using Base::begin;
    using Base::rbegin;
    using Base::end;
    using Base::rend;
    using Base::clear;
    using Base::insert;
    using Base::find;
    using Base::empty;
    using Base::count;

    using typename Base::iterator;
    using typename Base::const_iterator;

    using typename Base::mapped_type;
    using typename Base::value_type;
    //@}

    ///Test whether the map contains the given key.
    inline bool has(const Key& key) const
    {
      return Base::find(key) != Base::end();
    }

    /**
      @brief Return a constant reference to the element whose key is @p key.

      @exception IllegalKey if the given key does not exist
    */
    const T& operator[](const Key& key) const;

    /// Return a mutable reference to the element whose key is @p key. If an element with the key @p key does not exist, it is inserted.
    T& operator[](const Key& key);


    inline bool equals(const Map<Key, T>& other) const
    {
      return operator==(
        static_cast<typename Map<Key, T>::Base>(*this),
        static_cast<typename Map<Key, T>::Base>(other)
        );
    }

  };

  //******************************************************************************************
  // Implementations of template methods
  //******************************************************************************************

  template <class Key, class T>
  const T& Map<Key, T>::operator[](const Key& key) const
  {
    ConstIterator it = this->find(key);
    if (it == Base::end())
    {
      throw IllegalKey(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    else
    {
      return it->second;
    }
  }

  template <class Key, class T>
  T& Map<Key, T>::operator[](const Key& key)
  {
    Iterator it = this->find(key);
    if (it == Base::end())
    {
      it = this->insert(ValueType(key, T())).first;
    }
    return it->second;
  }

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_MAP_H
