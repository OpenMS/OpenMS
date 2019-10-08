// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  /** @brief A base class defining a common interface for all classes having a unique id.

   Have a look at RichPeak2D for an example how to extend a class to support unique ids.

   @sa UniqueIdGenerator, RichPeak2D

   @ingroup Concept
   */
  class OPENMS_DLLAPI UniqueIdInterface
  {
public:

    /**@brief This is the invalid unique id (cast it to a UInt64 if you like)

     It is represented as an \c enum because \c static class members lead to bugs and linker errors all the time...
     */
    enum
    {
      INVALID = 0
    };

    /**@brief Returns true if the unique_id is valid, false otherwise.

     Currently, an invalid unique id is represented by UInt64(0), but please prefer using this method for clarity.
     */
    inline static bool
    isValid(UInt64 unique_id)
    {
      return unique_id != INVALID;
    }

    /// Default constructor - the unique id will be <i>invalid</i>
    UniqueIdInterface() :
      unique_id_(UInt64(INVALID))
    {
      //do not use clearUniqueId(), as it will be slower and creates a warranted valgrind warning;
    }

    /// Copy constructor - copies the unique id
    UniqueIdInterface(const UniqueIdInterface & rhs) = default;

    /// Move constructor 
    UniqueIdInterface(UniqueIdInterface && rhs) = default;

    /// Assignment operator - copies the unique id
    UniqueIdInterface & operator=(UniqueIdInterface const & rhs) = default;

    /// Move Assignment operator - copies the unique id
    UniqueIdInterface& operator=(UniqueIdInterface&&) & = default;

    /// Destructor
    virtual ~UniqueIdInterface() = default;

    /// Equality comparison operator - the unique ids must be equal (!)
    bool
    operator==(UniqueIdInterface const & rhs) const
    {
      return unique_id_ == rhs.unique_id_;
    }

    /// Non-mutable access to unique id - returns the unique id.
    UInt64
    getUniqueId() const
    {
      return unique_id_;
    }

    /// Clear the unique id.  The new unique id will be invalid.  Returns 1 if the unique id was changed, 0 otherwise.
    Size
    clearUniqueId()
    {
      if (hasValidUniqueId())
      {
        unique_id_ = 0;
        return 1;
      }
      else
        return 0;
    }

    void
    swap(UniqueIdInterface & from)
    {
      std::swap(unique_id_, from.unique_id_);
      return;
    }

    /// Returns whether the unique id is valid.  Returns 1 if the unique id is valid, 0 otherwise.
    Size
    hasValidUniqueId() const
    {
      return isValid(unique_id_);
    }

    /// Returns whether the unique id is invalid.  Returns 1 if the unique id is invalid, 0 otherwise.
    Size
    hasInvalidUniqueId() const
    {
      return !isValid(unique_id_);
    }

    /// Assigns a new, valid unique id.  Always returns 1.
    Size
    setUniqueId()
    {
      unique_id_ = UniqueIdGenerator::getUniqueId();
      return 1;
    }

    /// Assigns a valid unique id, but only if the present one is invalid.  Returns 1 if the unique id was changed, 0 otherwise.
    Size
    ensureUniqueId()
    {
      if (!hasValidUniqueId())
      {
        unique_id_ = UniqueIdGenerator::getUniqueId();
        return 1;
      }
      else
        return 0;
    }

    /// Assigns the given unique id.
    void
    setUniqueId(UInt64 rhs)
    {
      unique_id_ = rhs;
    }

    /**@brief Mutable access to unique id.

     This is designed to work well with \c id attributes in XML, which are prefixed with letters.
     The portion of the String after the last underscore is extracted and parsed as a UInt64.  <i>It must consist of digits only.</i>
     For example, <code>some_feature.setUniqueId("f_12345_00067890")</code> is equivalent to <code>some_feature.setUniqueId(67890)</code>

     */
    void
    setUniqueId(const String & rhs);

protected:
    /// the unique id
    UInt64 unique_id_;

  };

} //namespace OpenMS

