// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>

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
    TextFile::load(filename, false, first_n, false, "#");
  }

  void CsvFile::load(const String& filename, char is, bool ie, Int first_n)
  {
    itemseperator_ = is;
    itemenclosed_ = ie;
    TextFile::load(filename, true, first_n, false, "#");
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

  bool CsvFile::getRow(Size row, StringList& list) const
  {
    // it is assumed that the value to be cast won't be so large to overflow an ulong int
    if (static_cast<int>(row) > static_cast<int>(TextFile::buffer_.size()) - 1)
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    bool splitted = buffer_[row].split(itemseperator_, list);
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
    return TextFile::buffer_.size();
  }

} // namespace OpenMS
