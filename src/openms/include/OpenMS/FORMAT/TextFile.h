// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

      @param filename @see load()
      @param trim_lines @see load()
      @param first_n @see load()
      @param skip_empty_lines @see load()

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    TextFile(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false);

    /**
      @brief Loads data from a text file.

      @param filename The input file name
      @param trim_lines Whether or not the lines are trimmed when reading them from file
      @param first_n If set, only @p first_n lines the lines from the beginning of the file are read
      @param skip_empty_lines Should empty lines be skipped? If used in conjunction with @p trim_lines, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false);

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
      @brief Platform-agnostic getline() which can deal with all line endings (\r, \r\n, \n)

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

