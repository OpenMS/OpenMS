// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

namespace OpenMS
{
  /**
      @brief This class provides some basic file handling methods for text files.

  @ingroup FileIO
  */
  class OPENMS_DLLAPI TextFile
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


    ///Default constructor
    TextFile();

    /// destructor
    virtual ~TextFile();

    /**
      @brief Constructor with filename

      @param filename The input file name
      @param trim_lines Whether or not the lines are trimmed when reading them from file
      @param first_n If set, only @p first_n lines the lines from the beginning of the file are read
      @param skip_empty_lines Should empty lines be skipped? If used in conjunction with @p trim_lines, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines.
      @param comment_symbol Lines prefixed with this string are skipped. Comment lines do not count towards the total number of read lines.
      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    TextFile(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false, const String& comment_symbol = "");

    /**
      @brief Loads data from a text file into the internal buffer.

      Retrieve the data using begin() and end().

      @param filename The input file name
      @param trim_lines Whether or not the lines are trimmed when reading them from file
      @param first_n If set, only @p first_n lines the lines from the beginning of the file are read
      @param skip_empty_lines Should empty lines be skipped? If used in conjunction with @p trim_lines, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines.
      @param comment_symbol Lines prefixed with this string are skipped. Comment lines do not count towards the total number of read lines.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false, const String& comment_symbol = "");

    /**
      @brief Writes the data to a file

      @note This function uses platform-dependent line breaks

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename);

    /// Operator for appending entries with less code
    template <typename StringType>
    TextFile& operator<<(const StringType& string)
    {
      buffer_.push_back(static_cast<String>(string));
      return *this;
    }

    template <typename StringType>
    void addLine(const StringType& line)
    {
      buffer_.push_back(static_cast<String>(line));
    }

    /**
      @brief Platform-agnostic getline() which can deal with all line endings (\\r, \\r\\n, \\n)

      Line endings will be removed from the resulting string.
    
    */
    static std::istream& getLine(std::istream& is, std::string& t);

    /**
      @brief Gives access to the underlying text buffer.
    */
    ConstIterator begin() const;

    Iterator begin();
    /**
     @brief Gives access to the underlying text buffer.
     */
    ConstIterator end() const;

    Iterator end();

protected:
    /// Internal buffer storing the lines before writing them to the file.
    std::vector<String> buffer_;
  };

} // namespace OpenMS

