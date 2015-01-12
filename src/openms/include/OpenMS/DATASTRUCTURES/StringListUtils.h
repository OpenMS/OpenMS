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

#ifndef OPENMS_DATASTRUCTURES_STRINGLIST_H
#define OPENMS_DATASTRUCTURES_STRINGLIST_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/OpenMSConfig.h>

class QStringList;

namespace OpenMS
{

  /**
    @brief Utilities operating on lists of Strings

    @ingroup Datastructures
  */
  class OPENMS_DLLAPI StringListUtils
  {
public:
    /** @name Type definitions
    */
    //@{
    /// Mutable iterator
    typedef std::vector<String>::iterator Iterator;
    /// Non-mutable iterator
    typedef std::vector<String>::const_iterator ConstIterator;
    /// Mutable reverse iterator
    typedef std::vector<String>::reverse_iterator ReverseIterator;
    /// Non-mutable reverse iterator
    typedef std::vector<String>::const_reverse_iterator ConstReverseIterator;
    //@}

    /// Creates a StringList from a QStringList
    static StringList fromQStringList(const QStringList& rhs);

    ///@name Search methods
    //@{
    /**
      @brief Searches for the first line that starts with @p text beginning at line @p start

      @param start Iterator pointing to the initial position to search. (note: that this does not need to correspond to the beginning of the container.
      @param end Iterator pointing to the end final position of the sequence to search.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static Iterator searchPrefix(const Iterator& start, const Iterator& end, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that starts with @p text beginning at line @p start

      @param start Iterator pointing to the initial position to search. (note: that this does not need to correspond to the beginning of the container.
      @param end Iterator pointing to the end final position of the sequence to search.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static ConstIterator searchPrefix(const ConstIterator& start, const ConstIterator& end, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that starts with @p text in the StringList @p container.

      @param container The StringList that should be searched.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static ConstIterator searchPrefix(const StringList& container, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that starts with @p text in the StringList @p container.

      @param container The StringList that should be searched.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static Iterator searchPrefix(StringList& container, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that ends with @p text beginning at line @p start

      @param start Iterator pointing to the initial position to search. (note: that this does not need to correspond to the beginning of the container.
      @param end Iterator pointing to the end final position of the sequence to search.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static Iterator searchSuffix(const Iterator& start, const Iterator& end, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that ends with @p text beginning at line @p start

      @param start Iterator pointing to the initial position to search. (note: that this does not need to correspond to the beginning of the container.
      @param end Iterator pointing to the end final position of the sequence to search.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static ConstIterator searchSuffix(const ConstIterator& start, const ConstIterator& end, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that ends with @p text in the StringList @p container.

      @param container The StringList that should be searched.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static ConstIterator searchSuffix(const StringList& container, const String& text, bool trim = false);

    /**
      @brief Searches for the first line that ends with @p text in the StringList @p container.

      @param container The StringList that should be searched.
      @param text The text to find
      @param trim Whether the line is trimmed before
      @return Returns an iterator to the matching entry. If no line matches end is returned.
    */
    static Iterator searchSuffix(StringList& container, const String& text, bool trim = false);


    //@}

    /**
      @brief Transforms all strings contained in the passed StringList to upper case.

      @param sl The StringList to convert to upper case.
    */
    static void toUpper(StringList& sl);

    /**
      @brief Transforms all strings contained in the passed StringList to lower case.

      @param The StringList to convert to lower case.
    */
    static void toLower(StringList& sl);

private:
    /// @cond INTERNAL
    struct TrimmableStringPredicate_
    {
      TrimmableStringPredicate_(const String& target, const bool trim) :
        trim_(trim),
        target_(target)
      {
        if (trim_) target_.trim();
      }

      inline String getValue(const String& value)
      {
        if (trim_)
        {
          // trim is not a const function so we need to create a copy first
          String cp = value;
          return cp.trim();
        }
        else
        {
          return value;
        }
      }

protected:
      /// Should the strings be trimmed.
      bool trim_;
      /// The target value that should be found.
      String target_;
    };

    /// Predicate to search in a StringList for a specific prefix.
    struct PrefixPredicate_ :
      TrimmableStringPredicate_
    {
      PrefixPredicate_(const String& target, const bool trim) :
        TrimmableStringPredicate_(target, trim)
      {}

      /**
        @brief Returns true if the (trimmed) value has the prefix @p target_.

        @param value The value to test.
        @return true if value has prefix target, false otherwise.
      */
      inline bool operator()(const String& value)
      {
        return getValue(value).hasPrefix(target_);
      }

    };

    /// Predicate to search in a StringList for a specific suffix.
    struct SuffixPredicate_ :
      TrimmableStringPredicate_
    {
      SuffixPredicate_(const String& target, const bool trim) :
        TrimmableStringPredicate_(target, trim)
      {}

      /**
       @brief Returns true if the (trimmed) value has the suffix @p target_.

       @param value The value to test.
       @return true if value has suffix target, false otherwise.
       */
      inline bool operator()(const String& value)
      {
        return getValue(value).hasSuffix(target_);
      }

    };
    /// @endcond INTERNAL

    /// hide c'tors to avoid instantiation of utils class
    StringListUtils() {}
    StringListUtils(const StringListUtils&){}
    StringListUtils& operator=(StringListUtils&){return *this; }
  };

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_STRINGLIST_H
