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

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{

  CsvFile::CsvFile() :
    TextFile(), itemseperator_(','), itemenclosed_(false)
  {

  }

  CsvFile::~CsvFile()
  {
  }

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
    // it is assumed that the value to be casted won't be so large to overflow an int
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
