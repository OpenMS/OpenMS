// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( push )
#pragma warning( disable : 4251 )     // disable MSVC dll-interface warning
#endif

class QStringList;

namespace OpenMS
{
  /**
      @brief String list

      This class is based on std::vector<String> but adds some methods for convenience.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI StringList :
    public std::vector<String>
  {
public:

    /** @name Type definitions
    */
    //@{
    /// Mutable iterator
    typedef iterator    Iterator;
    /// Non-mutable iterator
    typedef const_iterator  ConstIterator;
    /// Mutable reverse iterator
    typedef reverse_iterator    ReverseIterator;
    /// Non-mutable reverse iterator
    typedef const_reverse_iterator  ConstReverseIterator;
    //@}

    ///@name Constructors and assignment operators
    //@{
    /// Default constructor
    StringList();
    /// Copy constructor
    StringList(const StringList & rhs);
    /// Constructor from vector<String>
    StringList(const std::vector<String> & rhs);
    /// Constructor from vector<string>
    StringList(const std::vector<std::string> & rhs);
    /// Constructor from QStringList
    StringList(const QStringList & rhs);
    ///  Assignment operator
    StringList & operator=(const StringList & rhs);
    ///  Assignment operator from vector<String>
    StringList & operator=(const std::vector<String> & rhs);
    ///  Assignment operator vector<string>
    StringList & operator=(const std::vector<std::string> & rhs);
    //@}

    ///@name Search methods
    //@{
    /**
@brief Searches for the first line that starts with @p text beginning at line @p start

@param start the line to start the search in
@param text the text to find
@param trim whether the line is trimmed before
@return returns an iterator to the matching line. If no line matches, end() is returned
*/
    Iterator search(const Iterator & start, const String & text, bool trim = false);

    /**
        @brief Searches for the first line that starts with @p text

        This is an overloaded member function, provided for convenience.<br>
        It behaves essentially like the above function but the search is start at the beginning of the file
*/
    Iterator search(const String & text, bool trim = false);

    /**
    @brief Searches for the first line that ends with @p text beginning at line @p start

    @param start the line to start the search in
    @param text the text to find
    @param trim whether the line is trimmed before
    @return returns an iterator to the matching line. If no line matches, end() is returned
*/
    Iterator searchSuffix(const Iterator & start, const String & text, bool trim = false);

    /**
        @brief Searches for the first line that ends with @p text

        This is an overloaded member function, provided for convenience.

        It behaves essentially like searchSuffix(const Iterator&, const String&, bool) but the search starts at the beginning of the file
*/
    Iterator searchSuffix(const String & text, bool trim = false);

    /**
      @brief Searches for the first line that starts with @p text beginning at line @p start

      @param start the line to start the search in
      @param text the text to find
      @param trim whether the line is trimmed before
      @return returns an iterator to the matching line. If no line matches, end() is returned
    */
    ConstIterator search(const ConstIterator & start, const String & text, bool trim = false) const;

    /**
      @brief Searches for the first line that starts with @p text

      This is an overloaded member function, provided for convenience.<br>
      It behaves essentially like the above function but the search is start at the beginning of the file
    */
    ConstIterator search(const String & text, bool trim = false) const;

    /**
      @brief Searches for the first line that ends with @p text beginning at line @p start

      @param start the line to start the search in
      @param text the text to find
      @param trim whether the line is trimmed before
      @return returns an iterator to the matching line. If no line matches, end() is returned
    */
    ConstIterator searchSuffix(const ConstIterator & start, const String & text, bool trim = false) const;

    /**
      @brief Searches for the first line that ends with @p text

      This is an overloaded member function, provided for convenience.

      It behaves essentially like searchSuffix(const Iterator&, const String&, bool) but the search starts at the beginning of the file
    */
    ConstIterator searchSuffix(const String & text, bool trim = false) const;
    //@}

    ///Operator for appending entries with less code
    template <typename StringType>
    StringList & operator<<(const StringType & string)
    {
      this->push_back(string);
      return *this;
    }

    /// Returns a list that is created by splitting the given (comma-separated) string (String are not trimmed!)
    static StringList create(const String & list, const char splitter = ',');
    /// Returns a list that is created from an array of char*
    static StringList create(const char * const * list, UInt size);
    /// Returns if a string is contained in the list
    bool contains(const String & s) const;
    /// Transforms all contained strings to upper case
    void toUpper();
    /// Transforms all contained strings to lower case
    void toLower();

    /// Concatenate the string elements and putting the @p glue string between elements
    String concatenate(const String & glue = "") const;

    /// output stream operator
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const StringList & p);

  };

} // namespace OPENMS

#ifdef OPENMS_COMPILER_MSVC
#pragma warning( pop )
#endif

#endif // OPENMS_DATASTRUCTURES_STRINGLIST_H
