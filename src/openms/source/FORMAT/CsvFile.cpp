// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace std;

namespace OpenMS
{

  CsvFile::CsvFile() :
    TextFile(), itemseperator_(','), itemenclosed_(false)
  {

  }

  CsvFile::~CsvFile() = default;

  CsvFile::CsvFile(const String& filename, char is, bool ie, Int first_n) :
    TextFile(), itemseperator_(is), itemenclosed_(ie)
  {
    TextFile::load(filename, false, first_n);
  }

  void CsvFile::load(const String& filename, char is, bool ie, Int first_n)
  {
    itemseperator_ = is;
    itemenclosed_ = ie;
    TextFile::load(filename, true, first_n);
  }

  void CsvFile::store(const String& filename)
  {
    TextFile::store(filename);
  }

  void CsvFile::addRow(const StringList& list)
  {
    StringList elements = list;
    if (itemenclosed_)
    {
      for (Size i = 0; i < elements.size(); ++i)
      {
        elements[i].quote('"', String::NONE);
      }
    }
    String line;
    line.concatenate(elements.begin(), elements.end(), itemseperator_);
    addLine(line);
  }

  void CsvFile::clear()
  {
    buffer_.clear();
  }

  bool CsvFile::getRow(Size row, StringList& list)
  {
    // it is assumed that the value to be cast won't be so large to overflow an ulong int
    if (static_cast<int>(row) > static_cast<int>(CsvFile::rowCount()) - 1)
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // look for lines starting with '#' and increase the position of the retrieved row
    // this way 'content_row' will ignore comment rows for the count
    Size content_row = row;
    Size j = 0;
    // loop over all lines, looking for comments, until reaching the requested content 'row'
    while (j <= content_row)
    {
      // if the first letter is the comment char #, then add 1 to the retrieved row number
      if (buffer_[j].hasPrefix(Constants::COMMENT_CHAR))
      {
        ++content_row;
      }
      ++j;
    }

    bool splitted = buffer_[content_row].split(itemseperator_, list);
    if (!splitted)
    {
      return splitted;
    }
    for (Size i = 0; i < list.size(); i++)
    {
      if (itemenclosed_)
      {
        list[i] = list[i].substr(1, list[i].size() - 2);
      }
    }
    return true;
  }

  std::vector<String>::size_type CsvFile::rowCount() const
  {
    std::vector<String>::size_type content_size = 0;
    // loop over all lines and count, ignoring comments
    for (Size i = 0; i < TextFile::buffer_.size(); i++)
    {
      // if the first letter is not the comment char #, increase the count of content rows
      if (!buffer_[i].hasPrefix(Constants::COMMENT_CHAR))
      {
        ++content_size;
      }
    }
    return content_size;
  }

} // namespace OpenMS
