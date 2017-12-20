// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CSVFILE_H
#define OPENMS_FORMAT_CSVFILE_H

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

namespace OpenMS
{
  /**
    @brief This class handles csv files. Currently only loading is implemented.

    @note items are allowed to be enclosed by only one character e.g. "item" where " is enclosing character

    @ingroup FileIO
  */
  class OPENMS_DLLAPI CsvFile :
    private TextFile
  {
public:

    ///Default constructor
    CsvFile();

    /// destructor
    ~CsvFile() override;

    /**
      @brief Constructor with filename

      @param  filename The input file name.
      @param  is character which separates the items.
      @param  ie Whether or not every item is enclosed.
      @param  first_n Only the given number of lines are read, starting from the beginning of the file.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    CsvFile(const String& filename, char is = ',', bool ie = false, Int first_n = -1);

    /**
      @brief Loads data from a text file.

      @param  filename The input file name.
      @param  is character which separates the items.
      @param  ie Whether or not every item is enclosed.
      @param  first_n Only the given number of lines are read, starting from the beginning of the file.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const String& filename, char is = ',', bool ie = false, Int first_n = -1);

    /**
      @brief Stores the buffer's content into a file.

      @param filename The output filename.
    */
    void store(const String& filename);

    /**
      @brief Add a row to the buffer.

      @param list StringList which will contain all items of the row to add
    */
    void addRow(const StringList& list);

    /**
      @brief Clears the buffer

      Clears TextFile::buffer_
    */
    void clear();

    /**
      @brief writes all items from a row to list

      @param row the row which will be read
      @param list StringList which will contain all items of the row

      @exception Exception::InvalidIterator is thrown if the row is not existing

      @return  returns false if the given row could not be separated into items
    */
    bool getRow(Size row, StringList& list);

    /**
      @brief Returns the number of rows that were loaded from the file.

      @return The number of loaded rows.
    */
    std::vector<String>::size_type rowCount() const;

private:
    char itemseperator_;
    bool itemenclosed_;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_CSVFILE_H
