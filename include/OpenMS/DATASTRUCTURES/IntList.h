// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_INTLIST_H
#define OPENMS_DATASTRUCTURES_INTLIST_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( push )
#pragma warning( disable : 4251 )     // disable MSVC dll-interface warning
#endif

namespace OpenMS
{
  /**
      @brief Int list

      This class is based on std::vector<Int> but adds some methods for convenience.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI IntList :
    public std::vector<Int>
  {
public:

    ///@name Constructors and assignment operators
    //@{
    /// Default constructor
    IntList();
    /// Copy constructor
    IntList(const IntList & rhs);
    /// Constructor from vector<UInt>
    IntList(const std::vector<UInt> & rhs);
    /// Constructor from vector<Int>
    IntList(const std::vector<Int> & rhs);
    ///  Assignment operator
    IntList & operator=(const IntList & rhs);
    ///  Assignment operator from vector<Int>
    IntList & operator=(const std::vector<Int> & rhs);
    ///  Assignment operator from vector<UInt>
    IntList & operator=(const std::vector<UInt> & rhs);
    //@}

    ///Operator for appending entries with less code
    template <typename IntType>
    IntList & operator<<(IntType value)
    {
      this->push_back(value);
      return *this;
    }

    /// Returns a list that is created by splitting the given comma-separated string (String are not trimmed!)
    static IntList create(const String & list);
    ///Returns a list that is created by converting every string element of the given StringList
    static IntList create(const StringList & list);
    /// Returns if a string is contains in the list
    bool contains(Int s) const;


    /// output stream operator
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const IntList & p);

  };


} // namespace OPENMS

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( pop )
#endif

#endif // OPENMS_DATASTRUCTURES_INTLIST_H
