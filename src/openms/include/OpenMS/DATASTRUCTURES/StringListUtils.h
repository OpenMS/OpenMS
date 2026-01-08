// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/OpenMSConfig.h>

#include <QtCore/qcontainerfwd.h> // for QStringList

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

      @param sl The StringList to convert to lower case.
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

      inline String getValue(const String& value) const
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

