// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

namespace OpenMS
{
  /**
    @brief This class handles csv files. Currently only loading is implemented. Does NOT support comment lines!

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
    bool getRow(Size row, StringList& list) const;

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

